/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "functionObjectSurface.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkCellData.h"
#include "vtkCellDataToPointData.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositePolyDataMapper.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"

// VTK Readers
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeName(functionObjectSurface);
    addToRunTimeSelectionTable(surface, functionObjectSurface, dictionary);
}
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

static vtkSmartPointer<vtkPolyData> getPolyDataFile(const Foam::fileName& fName)
{
    // Not extremely elegant...
    vtkSmartPointer<vtkPolyData> dataset;

    if ("vtk" == fName.ext())
    {
        auto reader = vtkSmartPointer<vtkPolyDataReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }

    if ("vtp" == fName.ext())
    {
        auto reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

        reader->SetFileName(fName.c_str());
        reader->Update();
        dataset = reader->GetOutput();

        return dataset;
    }

    return dataset;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::functionObjectSurface::
functionObjectSurface
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    geometrySurface(parent, dict, colours, List<fileName>()),
    functionObjectBase(parent, dict, colours)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::functionObjectSurface::
addGeometry
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return false;
    }

    DebugInfo
        << "    Find surface " << functionObjectName_ << endl;

    const polySurface* surf =
    (
        geometryBase::parent_.storedObjects()
       .cfindObject<polySurface>(functionObjectName_)
    );

    // Treat surface with no faces/points like a missing surface
    surf = ((surf && surf->nPoints()) ? surf : nullptr);

    bool hasSurface = surf;


    // Retrieve the field association (CELL, POINT) for the given field

    unsigned fieldAssociation(0u);
    if (surf)
    {
        unsigned queried = surf->queryFieldAssociation(fieldName_);

        if (queried & polySurface::FACE_DATA)
        {
            fieldAssociation |= FieldAssociation::CELL_DATA;
        }
        if (queried & polySurface::POINT_DATA)
        {
            fieldAssociation |= FieldAssociation::POINT_DATA;
        }
    }

    // Reduce the information
    if (Pstream::parRun())
    {
        if (!hasSurface)
        {
            // No geometry - set all field association bits ON to ensure
            // it does not affect bitwise reduction.
            fieldAssociation = (~0u);
        }

        reduce(hasSurface, orOp<bool>());
        reduce(fieldAssociation, bitAndOp<unsigned>());
    }

    if (!hasSurface)
    {
        WarningInFunction
            << "No functionObject surface, or has no faces: "
            << functionObjectName_
            << endl;

        DebugInfo
            << "    Available surfaces:" << nl
            << geometryBase::parent_.storedObjects()
                .sortedNames<polySurface>() << endl;

        return false;
    }

    //// Pout<< "local surface = " << (surf ? surf->nFaces() : 0) << nl;


    // Create a vtkMultiPieceDataSet with vtkPolyData on the leaves
    vtkSmartPointer<vtkMultiPieceDataSet> multiPiece;

    // Requesting glyphs on the surface AND only have face data?
    // - just use the faceCentres directly and attach fields as CellData
    //   (not PointData).

    if
    (
        representation_ == rtGlyph
     && (fieldAssociation == FieldAssociation::CELL_DATA)
    )
    {
        multiPiece = gatherFaceCentres(surf);
    }
    else
    {
        multiPiece = gatherSurfacePieces(surf);
    }


    // Add the field (the information is consistent after last reduction).

    // Need field(s) for glyphs or colourByField:

    if (representation_ == rtGlyph || colourBy_ == cbField)
    {
        if (fieldAssociation == FieldAssociation::CELL_DATA)
        {
            addDimField<polySurfaceGeoMesh>
            (
                multiPiece,
                surf,
                fieldName_
            );
        }
        else if (fieldAssociation & FieldAssociation::POINT_DATA)
        {
            addDimField<polySurfacePointGeoMesh>
            (
                multiPiece,
                surf,
                fieldName_
            );
        }
    }


    // Now have a multi-piece dataset that is one of the following:
    //
    // - one-piece per processor (OpenFOAM = parallel, VTK=parallel)
    // - all pieces on master only (OpenFOAM = parallel, VTK=serial)

    // Re-query field information - we may have stored it differently
    // than the original source.

    fieldSummary fieldInfo = queryFieldSummary(fieldName_, multiPiece);
    fieldInfo.reduce();

    DebugInfo
        << "    Field " << fieldName_ << ' ' << fieldInfo.info() << endl;


    // Not rendered on this processor?
    // This is where we stop, but could also have an MPI barrier
    if (!renderer)
    {
        return true;
    }


    // Rendering

    {
        auto polyData = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();

        polyData->SetInputData(multiPiece);
        polyData->Update();

        if (representation_ == rtGlyph)
        {
            addGlyphs
            (
                position,
                fieldName_, fieldInfo,  // scaling
                fieldName_, fieldInfo,  // colouring
                maxGlyphLength_,
                polyData->GetOutput(),
                surfaceActor_,
                renderer
            );
        }
        else
        {
            vtkSmartPointer<vtkCellDataToPointData> cellToPoint;

            // CellData - Need a cell->point filter
            if (smooth_ && !fieldInfo.hasPointData())
            {
                cellToPoint = vtkSmartPointer<vtkCellDataToPointData>::New();
                cellToPoint->SetInputData(multiPiece);

                polyData->SetInputConnection(cellToPoint->GetOutputPort());
            }
            else
            {
                polyData->SetInputData(multiPiece);
            }
            polyData->Update();


            if (!smooth_)
            {
                addFeatureEdges(renderer, polyData);
            }

            auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(polyData->GetOutputPort());

            setField
            (
                position,
                fieldName_,
                (
                    smooth_
                  ? FieldAssociation::POINT_DATA
                  : FieldAssociation(fieldInfo.association_)
                ),
                mapper,
                renderer
            );

            surfaceActor_->SetMapper(mapper);

            setRepresentation(surfaceActor_);

            renderer->AddActor(surfaceActor_);
        }
    }

    return true;
}


bool Foam::functionObjects::runTimePostPro::functionObjectSurface::
addGeometryFromFile
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return false;
    }

    vtkSmartPointer<vtkPolyData> polyData;

    bool good = true;

    // File reading is serial (master only)
    if (Pstream::master())
    {
        fileName fName = getFileName("file", fieldName_);

        if (fName.size())
        {
            polyData = getPolyDataFile(fName);

            if (!polyData || polyData->GetNumberOfPoints() == 0)
            {
                good = false;

                WarningInFunction
                    << "Could not read "<< fName << nl
                    << "Only VTK (.vtp, .vtk) files are supported"
                    << endl;
            }
            else
            {
                DebugInfo
                    << "    Resolved surface " << fName << endl;
            }
        }
        else
        {
            good = false;

            WarningInFunction
                << "Unable to read file name from function object "
                << functionObjectName_ << " for field " << fieldName_
                << ". Surface will not be processed"
                << endl;
        }
    }

    reduce(good, andOp<bool>());

    if (!good)
    {
        return false;
    }

    // Only render on master
    if (!renderer || !Pstream::master())
    {
        return true;
    }

    fieldSummary fieldInfo = queryFieldSummary(fieldName_, polyData);
    // No reduction (serial)

    DebugInfo
        << "    Field " << fieldName_ << ' ' << fieldInfo.info() << endl;


    // Render

    if (representation_ == rtGlyph)
    {
        addGlyphs
        (
            position,
            fieldName_, fieldInfo,  // scaling
            fieldName_, fieldInfo,  // colouring
            maxGlyphLength_,
            polyData,
            surfaceActor_,
            renderer
        );
    }
    else
    {
        addFeatureEdges(renderer, polyData);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polyData);

        setField
        (
            position,
            fieldName_,
            queryFieldAssociation(fieldName_, polyData),
            mapper,
            renderer
        );

        surfaceActor_->SetMapper(mapper);

        setRepresentation(surfaceActor_);

        renderer->AddActor(surfaceActor_);
    }

    return true;
}


void Foam::functionObjects::runTimePostPro::functionObjectSurface::
addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_)
    {
        return;
    }

    if (liveObject_)
    {
        // Live source
        if (addGeometry(position, renderer))
        {
            return;
        }

        WarningInFunction
            << "No functionObject live source, or is empty: "
            << functionObjectName_
            << " ... attempting with file source"
            << endl;
    }
    else
    {
        DebugInfo
            << "Using file source only" << nl;
    }

    // File source
    addGeometryFromFile(position, renderer);
}


bool Foam::functionObjects::runTimePostPro::functionObjectSurface::clear()
{
    if (functionObjectBase::clear())
    {
        // Even for a "live" data source we allow file cleanup
        // (eg, from a previous run, etc)
        return removeFile("file", fieldName_);
    }

    return false;
}


// ************************************************************************* //

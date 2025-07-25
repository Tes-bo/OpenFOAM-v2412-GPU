/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
#include "cuttingPlaneFilter.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkCellDataToPointData.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositePolyDataMapper.h"
#include "vtkCutter.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkPlane.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeName(cuttingPlaneFilter);
    addToRunTimeSelectionTable(surface, cuttingPlaneFilter, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::cuttingPlaneFilter::cuttingPlaneFilter
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    volumeFilter(parent, dict, colours),
    fieldVisualisationBase(dict, colours),
    plane_(dict),
    values_()
{
    dict.readIfPresent("offsets", values_);

    if (values_.empty())
    {
        values_.resize(1, scalar(0));
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::cuttingPlaneFilter::
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

    if (needsCollective())
    {
        Info<< type() << " : Not available for collective operation" << endl;
        return false;
    }

    DebugInfo << "    Adding cutting plane" << endl;


    // Bookkeeping for vtkUnstructuredGrid
    vtk::vtuAdaptor adaptor;
    vtkSmartPointer<vtkMultiPieceDataSet> multiPiece = mesh(adaptor);


    // Add (scalar/vector) field.
    // - Need field(s) for glyphs or colourByField:

    int nCmpt = 0;
    if (representation_ == rtGlyph || colourBy_ == cbField)
    {
        const auto* ioptr =
            parent().mesh().cfindObject<regIOobject>(fieldName_);

        if (!nCmpt)
        {
            nCmpt = addDimField<scalar>
            (
                multiPiece, adaptor, ioptr, fieldName_
            );
        }
        if (!nCmpt)
        {
            nCmpt = addDimField<vector>
            (
                multiPiece, adaptor, ioptr, fieldName_
            );
        }
    }


    // Now have a multi-piece dataset that is one of the following:
    //
    // - one-piece per processor (OpenFOAM = parallel, VTK=parallel)


    // Re-query field information - we may have stored it differently
    // than the original source.

    fieldSummary fieldInfo = queryFieldSummary(fieldName_, multiPiece);
    fieldInfo.reduce();


    // Not rendered on this processor?
    // This is where we stop, but could also have an MPI barrier
    if (!renderer)
    {
        return true;
    }


    // Rendering
    {
        // OpenFOAM plane -> vtkPlane definition

        auto pln = vtkSmartPointer<vtkPlane>::New();

        pln->SetNormal
        (
            plane_.normal().x(),
            plane_.normal().y(),
            plane_.normal().z()
        );
        pln->SetOrigin
        (
            plane_.origin().x(),
            plane_.origin().y(),
            plane_.origin().z()
        );


        // Plane cutting algorithm

        auto cutter = vtkSmartPointer<vtkCutter>::New();

        cutter->SetInputData(multiPiece);
        cutter->SetCutFunction(pln);

        cutter->SetNumberOfContours(values_.size());

        forAll(values_, pointi)
        {
            cutter->SetValue(pointi, values_[pointi]);
        }

        cutter->SetInputArrayToProcess
        (
            (nCmpt == 3 ? 1 : 0), // index: scalars(0), vectors(1)
            0, // port
            0, // connection
            vtkDataObject::FIELD_ASSOCIATION_CELLS,
            fieldName_.c_str()
        );

        cutter->Modified();
        cutter->Update();

        auto polyData = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();

        polyData->SetInputConnection(cutter->GetOutputPort());
        polyData->Update();

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(polyData->GetOutputPort());

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
                cellToPoint->SetInputConnection(cutter->GetOutputPort());

                polyData->SetInputConnection(cellToPoint->GetOutputPort());
                polyData->Update();
            }

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


void Foam::functionObjects::runTimePostPro::cuttingPlaneFilter::
addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (visible_)
    {
        // Live source
        if (addGeometry(position, renderer))
        {
            return;
        }

        WarningInFunction
            << "Unsupported for OpenFOAM parallel and VTK serial"
            << endl;
    }
}


bool Foam::functionObjects::runTimePostPro::cuttingPlaneFilter::clear()
{
    return true;
}


// ************************************************************************* //

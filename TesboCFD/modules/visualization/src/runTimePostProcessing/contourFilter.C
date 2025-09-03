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
#include "contourFilter.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkCellDataToPointData.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositePolyDataMapper.h"
#include "vtkContourFilter.h"
#include "vtkMultiPieceDataSet.h"
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
    defineTypeName(contourFilter);
    addToRunTimeSelectionTable(surface, contourFilter, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::contourFilter::contourFilter
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    volumeFilter(parent, dict, colours),
    fieldVisualisationBase(dict, colours),
    colourFieldName_(dict.get<word>("colourField")),
    values_()
{
    dict.readEntry("values", values_);

    // Extra safety
    if (values_.empty())
    {
        values_.resize(1, scalar(0));
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::contourFilter::addGeometry
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

    DebugInfo << "    Adding iso-surface" << endl;

    // Bookkeeping for vtkUnstructuredGrid
    vtk::vtuAdaptor adaptor;
    vtkSmartPointer<vtkMultiPieceDataSet> multiPiece = mesh(adaptor);


    // Add (scalar/vector) field.
    // - always need field(s) for glyphs or colourByField:

    int nCmpt = 0;
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


    // If the input is vector, need magnitude

    word magFieldName = fieldName_;

    if (nCmpt == 3)
    {
        addMagField(fieldName_, multiPiece);
        magFieldName = "mag(" + fieldName_ + ")";
    }

    // Colouring
    nCmpt = 0;
    if (colourBy_ == cbField && fieldName_ != colourFieldName_)
    {
        const auto* ioptr =
            parent().mesh().cfindObject<regIOobject>(fieldName_);

        if (!nCmpt)
        {
            nCmpt = addDimField<scalar>
            (
                multiPiece, adaptor, ioptr, colourFieldName_
            );
        }
        if (!nCmpt)
        {
            nCmpt = addDimField<vector>
            (
                multiPiece, adaptor, ioptr, colourFieldName_
            );
        }
    }


    // Now have a multi-piece dataset that is one of the following:
    //
    // - one-piece per processor (OpenFOAM = parallel, VTK=parallel)


    // Re-query field information - we may have stored it differently
    // than the original source.

    fieldSummary fieldInfo = queryFieldSummary(magFieldName, multiPiece);
    fieldInfo.reduce();

    fieldSummary colourFieldInfo =
        queryFieldSummary(colourFieldName_, multiPiece);
    colourFieldInfo.reduce();


    DebugInfo
        << "    Field " << fieldName_ << ' ' << fieldInfo.info() << nl
        << "    Field " << colourFieldName_ << ' ' << colourFieldInfo.info()
        << endl;


    // Not rendered on this processor?
    // This is where we stop, but could also have an MPI barrier
    if (!renderer)
    {
        return true;
    }


    // Rendering
    {
        auto contour = vtkSmartPointer<vtkContourFilter>::New();

        vtkSmartPointer<vtkCellDataToPointData> cellToPoint;

        // CellData - Need a cell->point filter
        if (!fieldInfo.hasPointData() || !colourFieldInfo.hasPointData())
        {
            cellToPoint = vtkSmartPointer<vtkCellDataToPointData>::New();
            cellToPoint->SetInputData(multiPiece);

            contour->SetInputConnection(cellToPoint->GetOutputPort());
        }
        else
        {
            contour->SetInputData(multiPiece);
        }

        contour->SetNumberOfContours(values_.size());
        forAll(values_, valuei)
        {
            contour->SetValue(valuei, values_[valuei]);
        }

        contour->SetInputArrayToProcess
        (
            0,  // index: scalars(0)
            0,  // port
            0,  // connection
            vtkDataObject::FIELD_ASSOCIATION_POINTS,
            magFieldName.c_str()
        );

        contour->Modified();
        contour->Update();

        auto polyData = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();

        polyData->SetInputConnection(contour->GetOutputPort());
        polyData->Update();

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(polyData->GetOutputPort());

        if (representation_ == rtGlyph)
        {
            addGlyphs
            (
                position,
                colourFieldName_, colourFieldInfo,  // scaling
                colourFieldName_, colourFieldInfo,  // colouring
                maxGlyphLength_,
                polyData->GetOutput(),
                surfaceActor_,
                renderer
            );
        }
        else
        {
            setField
            (
                position,
                colourFieldName_,
                FieldAssociation::POINT_DATA,
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


void Foam::functionObjects::runTimePostPro::contourFilter::addGeometryToScene
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


bool Foam::functionObjects::runTimePostPro::contourFilter::clear()
{
    return true;
}


// ************************************************************************* //

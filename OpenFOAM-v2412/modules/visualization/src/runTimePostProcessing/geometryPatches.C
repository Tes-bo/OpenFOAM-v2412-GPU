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
#include "geometryPatches.H"
#include "fvMesh.H"
#include "volFields.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "foamVtkTools.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkCellArray.h"
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


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeName(geometryPatches);
    addToRunTimeSelectionTable(surface, geometryPatches, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometryPatches::geometryPatches
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    geometrySurface(parent, dict, colours, List<fileName>()),
    fieldVisualisationBase(dict, colours),
    selectPatches_(),
    nearCellValue_(dict.getOrDefault("nearCellValue", false))
{
    dict.readEntry("patches", selectPatches_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::geometryPatches::addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    if (!visible_ || selectPatches_.empty())
    {
        return;
    }

    const polyBoundaryMesh& pbm = parent().mesh().boundaryMesh();

    bitSet selectedPatchIds(pbm.size());

    for (const polyPatch& pp : pbm)
    {
        if (isType<emptyPolyPatch>(pp) || pp.empty())
        {
            continue;
        }
        else if (isA<processorPolyPatch>(pp))
        {
            break; // No processor patches
        }

        if (selectPatches_.match(pp.name()))
        {
            selectedPatchIds.set(pp.index());
        }
    }


    labelListList patchIds(Pstream::nProcs());
    patchIds[Pstream::myProcNo()] = selectedPatchIds.sortedToc();

    Pstream::allGatherList(patchIds);

    label nPieces = 0;
    for (const labelList& ids : patchIds)
    {
        nPieces += ids.size();
    }

    /// Pout<< "add patches: " << selectPatches_ << nl
    ///     << " = " << patchIds << " == " << nPieces << " total" << nl;

    if (!nPieces)
    {
        WarningInFunction
            << "No patches selected: " << flatOutput(selectPatches_)
            << endl;
        return;
    }

    DebugInfo << "    Add geometry patches" << nl;


    // Create a vtkMultiPieceDataSet with vtkPolyData on the leaves

    // When adding fields, only map scalar/vector fields.
    // - this is easier and all we mostly ever need

    vtkSmartPointer<vtkMultiPieceDataSet> multiPiece;

    // Requesting glyphs on the surface - just use the faceCentres directly
    // and attach fields as CellData (not PointData).

    if (representation_ == rtGlyph)
    {
        multiPiece = gatherPatchFaceCentres(patchIds);
    }
    else
    {
        multiPiece = gatherPatchPieces(patchIds);
    }


    // Add (scalar/vector) field.
    // - Need field(s) for glyphs or colourByField:

    int nCmpt = 0;
    if (representation_ == rtGlyph || colourBy_ == cbField)
    {
        if (!nCmpt)
        {
            nCmpt = addPatchField<scalar>
            (
                multiPiece,
                patchIds,
                parent().mesh().cfindObject<volScalarField>(fieldName_),
                fieldName_
            );
        }
        if (!nCmpt)
        {
            nCmpt = addPatchField<vector>
            (
                multiPiece,
                patchIds,
                parent().mesh().cfindObject<volVectorField>(fieldName_),
                fieldName_
            );
        }
    }

    // Now have a multi-piece dataset with
    // one piece per patch and processor.
    //
    // For VTK=parallel, these pieces reside on their original processors.
    // For VTK=serial, they are master only

    // Re-query actually field information, since we may have stored it
    // somewhere slightly different than the original source.

    fieldSummary fieldInfo = queryFieldSummary(fieldName_, multiPiece);
    fieldInfo.reduce();

    DebugInfo
        << "    Field " << fieldName_ << ' ' << fieldInfo.info() << endl;

    // Not rendering on this processor?
    // This is where we stop, but could also have a MPI barrier
    if (!renderer)
    {
        return;
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
                  : FieldAssociation::CELL_DATA
                ),
                mapper,
                renderer
            );

            surfaceActor_->SetMapper(mapper);

            setRepresentation(surfaceActor_);

            renderer->AddActor(surfaceActor_);
        }
    }
}


void Foam::functionObjects::runTimePostPro::geometryPatches::updateActors
(
    const scalar position
)
{
    if (!visible_)
    {
        return;
    }

    surface::updateActors(position);
}


bool Foam::functionObjects::runTimePostPro::geometryPatches::clear()
{
    return true;
}


// ************************************************************************* //

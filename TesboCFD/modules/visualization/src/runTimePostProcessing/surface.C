/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
#include "surface.H"
#include "runTimePostProcessing.H"

#include "foamVtkTools.H"
#include "polySurfaceFields.H"
#include "polySurfacePointFields.H"

// VTK includes
#include "vtkActor.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkFeatureEdges.h"
#include "vtkMultiPieceDataSet.h"
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
    defineTypeName(surface);
    defineRunTimeSelectionTable(surface, dictionary);
}
}
}


const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::surface::representationType
>
Foam::functionObjects::runTimePostPro::surface::representationTypeNames
({
    { representationType::rtNone, "none" },
    { representationType::rtGlyph, "glyph" },
    { representationType::rtWireframe, "wireframe" },
    { representationType::rtSurface, "surface" },
    { representationType::rtSurfaceWithEdges, "surfaceWithEdges" },
});


// * * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * //

// These need to shift elsewhere

vtkCellData* Foam::vtk::Tools::GetCellData(vtkDataSet* dataset)
{
    if (dataset) return dataset->GetCellData();
    return nullptr;
}

vtkPointData* Foam::vtk::Tools::GetPointData(vtkDataSet* dataset)
{
    if (dataset) return dataset->GetPointData();
    return nullptr;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::surface::setRepresentation
(
    vtkActor* actor
) const
{
    geometryBase::initialiseActor(actor);

    switch (representation_)
    {
        case rtNone:
        {
            actor->VisibilityOff();
            break;
        }
        case rtWireframe:
        {
            // Note: colour is set using general SetColor, not SetEdgeColor
            actor->GetProperty()->SetRepresentationToWireframe();
            break;
        }
        case rtGlyph:
        case rtSurface:
        {
            actor->GetProperty()->SetBackfaceCulling(backFaceCulling_);
            actor->GetProperty()->SetFrontfaceCulling(frontFaceCulling_);
            actor->GetProperty()->SetRepresentationToSurface();
            break;
        }
        case rtSurfaceWithEdges:
        {
            actor->GetProperty()->SetBackfaceCulling(backFaceCulling_);
            actor->GetProperty()->SetFrontfaceCulling(frontFaceCulling_);
            actor->GetProperty()->SetRepresentationToSurface();
            actor->GetProperty()->EdgeVisibilityOn();
            break;
        }
    }
}


void Foam::functionObjects::runTimePostPro::surface::addFeatureEdges
(
    vtkRenderer* renderer,
    vtkFeatureEdges* featureEdges
) const
{
    if (!featureEdges)
    {
        return;
    }

    featureEdges->BoundaryEdgesOn();
    featureEdges->FeatureEdgesOn();
    featureEdges->ManifoldEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    /// featureEdges->SetFeatureAngle(60);
    featureEdges->ColoringOff();
    featureEdges->Update();

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(featureEdges->GetOutputPort());
    mapper->ScalarVisibilityOff();

    edgeActor_->GetProperty()->SetSpecular(0);
    edgeActor_->GetProperty()->SetSpecularPower(20);
    edgeActor_->GetProperty()->SetRepresentationToWireframe();
    edgeActor_->SetMapper(mapper);

    renderer->AddActor(edgeActor_);
}


void Foam::functionObjects::runTimePostPro::surface::addFeatureEdges
(
    vtkRenderer* renderer,
    vtkPolyData* data
) const
{
    if (featureEdges_)
    {
        auto featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
        featureEdges->SetInputData(data);

        addFeatureEdges(renderer, featureEdges);
    }
}


void Foam::functionObjects::runTimePostPro::surface::addFeatureEdges
(
    vtkRenderer* renderer,
    vtkCompositeDataGeometryFilter* input
) const
{
    if (featureEdges_)
    {
        auto featureEdges = vtkSmartPointer<vtkFeatureEdges>::New();
        featureEdges->SetInputConnection(input->GetOutputPort());

        addFeatureEdges(renderer, featureEdges);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::surface::surface
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    geometryBase(parent, dict, colours),
    representation_
    (
        representationTypeNames.get("representation", dict)
    ),
    featureEdges_(dict.getOrDefault("featureEdges", false)),
    backFaceCulling_(dict.getOrDefault("backFaceCulling", false)),
    frontFaceCulling_(dict.getOrDefault("frontFaceCulling", false)),
    surfaceColour_(nullptr),
    edgeColour_(nullptr),
    surfaceActor_(),
    edgeActor_(),
    maxGlyphLength_(0)
{
    surfaceActor_ = vtkSmartPointer<vtkActor>::New();
    edgeActor_ = vtkSmartPointer<vtkActor>::New();

    surfaceColour_ = Function1<vector>::NewIfPresent("surfaceColour", dict);
    if (!surfaceColour_)
    {
        surfaceColour_.reset(colours["surface"]->clone().ptr());
    }

    edgeColour_ = Function1<vector>::NewIfPresent("edgeColour", dict);
    if (!edgeColour_)
    {
        edgeColour_.reset(colours["edge"]->clone().ptr());
    }

    if (representation_ == rtGlyph)
    {
        dict.readEntry("maxGlyphLength", maxGlyphLength_);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObjects::runTimePostPro::surface>
Foam::functionObjects::runTimePostPro::surface::New
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours,
    const word& surfaceType
)
{
    DebugInfo << "Selecting surface " << surfaceType << endl;

    auto* ctorPtr = dictionaryConstructorTable(surfaceType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "surface",
            surfaceType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<surface>(ctorPtr(parent, dict, colours));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::surface::~surface()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::surface::updateActors
(
    const scalar position
)
{
    if (!visible_)
    {
        return;
    }

    surfaceActor_->GetProperty()->SetOpacity(opacity(position));

    const vector sc = surfaceColour_->value(position);
    surfaceActor_->GetProperty()->SetColor(sc[0], sc[1], sc[2]);

    const vector ec = edgeColour_->value(position);
    surfaceActor_->GetProperty()->SetEdgeColor(ec[0], ec[1], ec[2]);

    if (featureEdges_)
    {
        vtkProperty* edgeProp = edgeActor_->GetProperty();

        edgeProp->SetLineWidth(2);
        edgeProp->SetOpacity(opacity(position));

        edgeProp->SetColor(ec[0], ec[1], ec[2]);
        edgeProp->SetEdgeColor(ec[0], ec[1], ec[2]);
    }
}


// ************************************************************************* //

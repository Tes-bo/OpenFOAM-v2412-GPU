/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
#include "geometryCloud.H"
#include "cloud.H"
#include "fvMesh.H"
#include "runTimePostProcessing.H"
#include "addToRunTimeSelectionTable.H"

// VTK includes
#include "vtkActor.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositePolyDataMapper.h"
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
    defineTypeName(geometryCloud);
    addToRunTimeSelectionTable(pointData, geometryCloud, dictionary);
}
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::geometryCloud::addCloudField
(
    vtkMultiPieceDataSet* multiPiece,
    const objectRegistry& obrTmp,
    const word& fieldName
) const
{
    const regIOobject* ioptr = obrTmp.cfindObject<regIOobject>(fieldName);

    return (multiPiece) &&
    (
        addCloudField<label>
        (
            multiPiece, ioptr, fieldName
        )
     || addCloudField<scalar>
        (
            multiPiece, ioptr, fieldName
        )
     || addCloudField<vector>
        (
            multiPiece, ioptr, fieldName
        )
     || addCloudField<sphericalTensor>
        (
            multiPiece, ioptr, fieldName
        )
     || addCloudField<symmTensor>
        (
            multiPiece, ioptr, fieldName
        )
     || addCloudField<tensor>
        (
            multiPiece, ioptr, fieldName
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometryCloud::geometryCloud
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    pointData(parent, dict, colours),
    fieldVisualisationBase(dict, colours),
    cloudName_(dict.get<word>("cloud")),
    colourFieldName_(dict.get<word>("colourField")),
    actor_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::geometryCloud::
~geometryCloud()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::runTimePostPro::geometryCloud::
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


    // This is similar (almost identical) to vtkCloud

    const auto* objPtr = parent().mesh().cfindObject<cloud>(cloudName_);
    if (!objPtr)
    {
        return false;
    }

    objectRegistry obrTmp
    (
        IOobject
        (
            "runTimePostPro::cloud::" + cloudName_,
            parent().mesh().time().constant(),
            parent().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    objPtr->writeObjects(obrTmp);

    const auto* pointsPtr = cloud::findIOPosition(obrTmp);

    if (!pointsPtr)
    {
        // This should be impossible
        return false;
    }


    // Create a vtkMultiPieceDataSet with vtkPolyData on the leaves
    auto multiPiece = gatherCloud(obrTmp);

    // Add in scaleField and colourField
    addCloudField(multiPiece, obrTmp, fieldName_);
    addCloudField(multiPiece, obrTmp, colourFieldName_);


    // Not rendered on this processor?
    // This is where we stop, but could also have an MPI barrier
    if (!renderer)
    {
        return true;
    }

    DebugInfo
        << "    Render cloud " << cloudName_ << endl;


    // Rendering

    actor_ = vtkSmartPointer<vtkActor>::New();

    {
        fieldSummary scaleFieldInfo =
            queryFieldSummary(fieldName_, multiPiece);

        fieldSummary colourFieldInfo =
            queryFieldSummary(colourFieldName_, multiPiece);

        auto polyData = vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();

        polyData->SetInputData(multiPiece);
        polyData->Update();

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        actor_->SetMapper(mapper);

        /// dataset->Print(std::cout);

        addGlyphs
        (
            position,
            fieldName_, scaleFieldInfo,         // scaling
            colourFieldName_, colourFieldInfo,  // colour
            maxGlyphLength_,
            polyData->GetOutput(),
            actor_,
            renderer
        );

        renderer->AddActor(actor_);
    }

    return true;
}


void Foam::functionObjects::runTimePostPro::geometryCloud::
addGeometryToScene
(
    const scalar position,
    vtkRenderer* renderer
)
{
    // Live source
    addGeometry(position, renderer);
}


void Foam::functionObjects::runTimePostPro::geometryCloud::updateActors
(
    const scalar position
)
{
    if (actor_)
    {
        const vector colour = pointColour_->value(position);

        vtkProperty* prop = actor_->GetProperty();

        prop->SetOpacity(opacity(position));

        prop->SetColor(colour[0], colour[1], colour[2]);
    }
}


bool Foam::functionObjects::runTimePostPro::geometryCloud::clear()
{
    return false;
}


// ************************************************************************* //

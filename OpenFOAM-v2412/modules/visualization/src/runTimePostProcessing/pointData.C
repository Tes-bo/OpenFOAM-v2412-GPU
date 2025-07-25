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
#include "pointData.H"
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkActor.h"
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
    defineTypeName(pointData);
    defineRunTimeSelectionTable(pointData, dictionary);
}
}
}

const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::pointData::representationType
>
Foam::functionObjects::runTimePostPro::pointData::representationTypeNames
({
    { representationType::rtSphere, "sphere" },
    { representationType::rtVector, "vector" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::pointData::addPoints
(
    const label framei,
    vtkActor* actor,
    vtkPolyDataMapper* mapper,
    vtkPolyData* data
) const
{
    geometryBase::initialiseActor(actor);

    vector colour = pointColour_->value(framei);
    actor->GetProperty()->SetColor(colour[0], colour[1], colour[2]);

    switch (representation_)
    {
        case rtSphere:
        case rtVector:
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::pointData::pointData
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
    maxGlyphLength_(dict.get<scalar>("maxGlyphLength")),
    pointColour_(nullptr)
{

    pointColour_ = Function1<vector>::NewIfPresent("pointColour", dict);
    if (!pointColour_)
    {
        pointColour_.reset(colours["point"]->clone().ptr());
    }

    switch (representation_)
    {
        case rtSphere:
        {
            break;
        }
        case rtVector:
        {
            break;
        }
    }

}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObjects::runTimePostPro::pointData>
Foam::functionObjects::runTimePostPro::pointData::New
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours,
    const word& pointDataType
)
{
    DebugInfo << "Selecting pointData " << pointDataType << endl;

    auto* ctorPtr = dictionaryConstructorTable(pointDataType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "pointData",
            pointDataType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<pointData>(ctorPtr(parent, dict, colours));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::pointData::~pointData()
{}


// ************************************************************************* //

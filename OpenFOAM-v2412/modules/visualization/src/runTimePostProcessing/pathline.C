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
#include "pathline.H"
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTubeFilter.h"
#include "vtkLookupTable.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{
    defineTypeName(pathline);
    defineRunTimeSelectionTable(pathline, dictionary);
}
}
}

const Foam::Enum
<
    Foam::functionObjects::runTimePostPro::pathline::representationType
>
Foam::functionObjects::runTimePostPro::pathline::representationTypeNames
({
    { representationType::rtNone, "none" },
    { representationType::rtLine, "line" },
    { representationType::rtTube, "tube" },
    { representationType::rtVector, "vector" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::runTimePostPro::pathline::addLines
(
    const label framei,
    vtkActor* actor,
    vtkPolyData* data
) const
{
    geometryBase::initialiseActor(actor);

    vector colour = lineColour_->value(framei);
    actor->GetProperty()->SetColor(colour[0], colour[1], colour[2]);

    auto* mapper = vtkPolyDataMapper::SafeDownCast(actor->GetMapper());

    switch (representation_)
    {
        case rtNone:
        {
            actor->VisibilityOff();
            break;
        }
        case rtLine:
        {
            mapper->SetInputData(data);
            mapper->Update();
            break;
        }
        case rtTube:
        {
            auto tubes = vtkSmartPointer<vtkTubeFilter>::New();
            tubes->SetInputData(data);
            tubes->SetRadius(tubeRadius_);
            tubes->SetNumberOfSides(20);
            tubes->CappingOn();
            tubes->Update();

            mapper->SetInputConnection(tubes->GetOutputPort());
            mapper->Update();
            break;
        }
        case rtVector:
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::pathline::pathline
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
    tubeRadius_(0.001),
    lineColour_(nullptr)
{
    lineColour_ = Function1<vector>::NewIfPresent("lineColour", dict);
    if (!lineColour_)
    {
        lineColour_.reset(colours["line"]->clone().ptr());
    }

    switch (representation_)
    {
        case rtNone:
        {
            break;
        }
        case rtLine:
        {
            break;
        }
        case rtTube:
        {
            dict.readEntry("tubeRadius", tubeRadius_);
            break;
        }
        case rtVector:
        {
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObjects::runTimePostPro::pathline>
Foam::functionObjects::runTimePostPro::pathline::New
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours,
    const word& pathlineType
)
{
    DebugInfo << "Selecting pathline " << pathlineType << endl;

    auto* ctorPtr = dictionaryConstructorTable(pathlineType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "pathline",
            pathlineType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<pathline>(ctorPtr(parent, dict, colours));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::pathline::~pathline()
{}


// ************************************************************************* //

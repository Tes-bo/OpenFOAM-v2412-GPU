/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "interfaceCompositionModel.H"
#include "phaseModel.H"
#include "phasePair.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseInter
{
    defineTypeNameAndDebug(interfaceCompositionModel, 0);
    defineRunTimeSelectionTable(interfaceCompositionModel, dictionary);
}
}

const Foam::Enum
<
    Foam::multiphaseInter::interfaceCompositionModel::modelVariable
>
Foam::multiphaseInter::interfaceCompositionModel::modelVariableNames_
{
    { modelVariable::T, "temperature" },
    { modelVariable::P, "pressure" },
    { modelVariable::Y, "massFraction" },
    { modelVariable::alpha, "alphaVolumeFraction" },
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseInter::interfaceCompositionModel::interfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    modelVariable_
    (
        modelVariableNames_.getOrDefault
        (
            "variable",
            dict,
            modelVariable::T
        )
    ),
    includeVolChange_(dict.getOrDefault("includeVolChange", true)),
    pair_(pair),
    speciesName_(dict.getOrDefault<word>("species", "none")),
    mesh_(pair_.from().mesh())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiphaseInter::interfaceCompositionModel>
Foam::multiphaseInter::interfaceCompositionModel::New
(
    const dictionary& dict,
    const phasePair& pair
)
{
    const word modelType
    (
        dict.get<word>("type")
      + "<"
      + pair.phase1().thermo().type()
      + ","
      + pair.phase2().thermo().type()
      + ">"
    );

    Info<< "Selecting interfaceCompositionModel for "
        << pair << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "interfaceCompositionModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, pair);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word Foam::multiphaseInter::interfaceCompositionModel
::transferSpecie() const
{
    return speciesName_;
}


const Foam::phasePair& Foam::multiphaseInter::interfaceCompositionModel
::pair() const
{
    return pair_;
}


const Foam::multiphaseInterSystem& Foam::multiphaseInter
::interfaceCompositionModel::fluid() const
{
    return pair().to().fluid();
}


const Foam::word& Foam::multiphaseInter::interfaceCompositionModel
::variable() const
{
    return modelVariableNames_[modelVariable_];
}


bool Foam::multiphaseInter::interfaceCompositionModel::includeDivU()
const noexcept
{
    return true;
}


bool Foam::multiphaseInter::interfaceCompositionModel::includeVolChange()
{
    return includeVolChange_;
}


// ************************************************************************* //

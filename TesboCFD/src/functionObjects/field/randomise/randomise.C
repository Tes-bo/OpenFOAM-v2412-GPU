/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "randomise.H"
#include "volFields.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(randomise, 0);
    addToRunTimeSelectionTable(functionObject, randomise, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::randomise::calcTemplate()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const auto* fieldPtr = cfindObject<VolFieldType>(fieldName_);

    if (fieldPtr)
    {
        const auto& field = *fieldPtr;

        resultName_ = scopedName(fieldName_ & "Random");

        auto trfield = tmp<VolFieldType>::New(field);
        auto& rfield = trfield.ref();

        Random rng(1234567);

        auto applyPerturbation = [&](Type& cellval)
        {
            Type rndPert;
            rng.randomise01(rndPert);
            rndPert = 2.0*rndPert - pTraits<Type>::one;
            rndPert /= mag(rndPert);

            cellval += magPerturbation_*rndPert;
        };

        if (this->volRegion::useAllCells())
        {
            for (Type& cellval : rfield)
            {
                applyPerturbation(cellval);
            }
        }
        else
        {
            for (const label celli : cellIDs())
            {
                applyPerturbation(rfield[celli]);
            }
        }

        return store(resultName_, trfield);
    }

    return false;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::randomise::calc()
{
    // Ensure volRegion is properly up-to-date.
    // Purge old fields if we need to etc.
    (void)volRegion::update();

    return
    (
        calcTemplate<scalar>()
     || calcTemplate<vector>()
     || calcTemplate<sphericalTensor>()
     || calcTemplate<symmTensor>()
     || calcTemplate<tensor>()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::randomise::randomise
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    volRegion(fieldExpression::mesh_, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::randomise::read(const dictionary& dict)
{
    fieldExpression::read(dict);
    volRegion::read(dict);

    dict.readEntry("magPerturbation", magPerturbation_);

    return true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "norm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(norm, 0);
    addToRunTimeSelectionTable(functionObject, norm, dictionary);
}
}


const Foam::Enum
<
    Foam::functionObjects::norm::normType
>
Foam::functionObjects::norm::normTypeNames
({
    { normType::L1 , "L1" },
    { normType::L2 , "L2" },
    { normType::LP , "Lp" },
    { normType::MAX , "max" },
    { normType::COMPOSITE , "composite" },
    { normType::FIELD , "divisorField" }
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Implementation
#include "normImpl.C"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::norm::calc()
{
    return
    (
        calcNorm<scalar>()
     || calcNorm<vector>()
     || calcNorm<sphericalTensor>()
     || calcNorm<symmTensor>()
     || calcNorm<tensor>()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::norm::norm
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    norm_(normType::L1),
    divisorPtr_(nullptr),
    divisorFieldName_(),
    p_(-1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::norm::read(const dictionary& dict)
{
    if (!fieldExpression::read(dict))
    {
        return false;
    }

    norm_ = normTypeNames.get("norm", dict);

    if (norm_ == normType::LP)
    {
        p_ = dict.getCheck<scalar>("p", scalarMinMax::ge(1));
    }

    if (norm_ == normType::COMPOSITE)
    {
        divisorPtr_ = Function1<scalar>::New("divisor", dict, &mesh_);

        if (!divisorPtr_)
        {
            FatalIOErrorInFunction(dict)
                << "The norm 'composite' needs the input entry 'divisor'."
                << abort(FatalIOError);
        }
    }

    if (norm_ == normType::FIELD)
    {
        divisorFieldName_ = dict.get<word>("divisorField");

        if (divisorFieldName_.empty())
        {
            FatalIOErrorInFunction(dict)
                << "The norm 'field' needs the input entry 'divisorField'."
                << abort(FatalIOError);
        }
    }

    return true;
}


// ************************************************************************* //

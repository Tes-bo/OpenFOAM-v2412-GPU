/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "multiply.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(multiply, 0);
    addToRunTimeSelectionTable(functionObject, multiply, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::multiply::multiply
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldsExpression(name, runTime, dict)
{
    setResultName("multiply");
}


bool Foam::functionObjects::multiply::calc()
{
    bool processed = false;

    Log << type() << ' ' << name() <<  " execute:" << nl;

    forAll(fieldNames_, i)
    {
        processed = false;

        if (i == 0)
        {
            initialiseResult<scalar>(fieldNames_[i]);
            initialiseResult<vector>(fieldNames_[i]);
            initialiseResult<sphericalTensor>(fieldNames_[i]);
            initialiseResult<symmTensor>(fieldNames_[i]);
            initialiseResult<tensor>(fieldNames_[i]);
        }
        else
        {
            multiplyResult<scalar>(fieldNames_[i], processed);
            multiplyResult<vector>(fieldNames_[i], processed);
            multiplyResult<sphericalTensor>(fieldNames_[i], processed);
            multiplyResult<symmTensor>(fieldNames_[i], processed);
            multiplyResult<tensor>(fieldNames_[i], processed);
        }
    }

    return processed;
}


// ************************************************************************* //

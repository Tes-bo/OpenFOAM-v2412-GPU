/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider
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

#include "exprResultStack.H"
#include "vector.H"
#include "tensor.H"
#include "symmTensor.H"
#include "sphericalTensor.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

    defineTypeNameAndDebug(exprResultStack, 0);
    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultStack,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        exprResult,
        exprResultStack,
        empty
    );

} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::exprResultStack::exprResultStack()
:
    expressions::exprResult()
{
    needsReset(true);  // Requires reset every timestep to work
}


Foam::expressions::exprResultStack::exprResultStack
(
    const exprResultStack& rhs
)
:
    expressions::exprResult(rhs)
{
    needsReset(true);  // Requires reset every timestep to work
}


Foam::expressions::exprResultStack::exprResultStack
(
    const dictionary& dict
)
:
    expressions::exprResult(dict)
{
    needsReset(true);  // Requires reset every timestep to work
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::expressions::exprResult
Foam::expressions::exprResultStack::pop()
{
    exprResult result;

    if (this->size() <= 0)
    {
        FatalErrorInFunction
            << "Trying to pop result from a empty queue" << endl
            << abort(FatalError);

        return result;
    }

    const bool ok =
    (
        popChecked<scalar>(result)
     || popChecked<vector>(result)
     || popChecked<tensor>(result)
     || popChecked<symmTensor>(result)
     || popChecked<sphericalTensor>(result)
    );

    if (!ok)
    {
        FatalErrorInFunction
            << "Unsupported value type " << valueType() << nl
            << abort(FatalError);
    }

    return result;
}


void Foam::expressions::exprResultStack::push(const exprResult& result)
{
    DebugInFunction << nl << "Pushing: " << result << nl;

    if (!hasValue())
    {
        // This is the first push
        exprResult::operator=(result);
    }
    else
    {
        if (valueType() != result.valueType())
        {
            FatalErrorInFunction
                << "Type of pushed value " << result.valueType()
                << " is not the expected type " << valueType() << nl
                << abort(FatalError);
        }

        const bool ok =
        (
            pushChecked<scalar>(result)
         || pushChecked<vector>(result)
         || pushChecked<tensor>(result)
         || pushChecked<symmTensor>(result)
         || pushChecked<sphericalTensor>(result)
        );

        if (!ok)
        {
            FatalErrorInFunction
                << "Unsupported value type " << valueType() << nl
                << abort(FatalError);
        }
    }

    DebugInFunction << "After push: " << *this << nl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::expressions::exprResultStack::operator=
(
    const exprResultStack& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    static_cast<exprResult&>(*this) = rhs;
}


void Foam::expressions::exprResultStack::operator=
(
    const exprResult& rhs
)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    DebugInFunction << nl;

    exprResult value
    (
        // Issue warning if the other result is not really uniform
        rhs.getUniform(1, false)
    );

    this->push(value);
}


// ************************************************************************* //

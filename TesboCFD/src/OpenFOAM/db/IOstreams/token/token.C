/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Note
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "token.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef token::compound tokenCompound;
    defineTypeName(tokenCompound);
    defineRunTimeSelectionTable(tokenCompound, empty);
}

const Foam::token Foam::token::undefinedToken;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::token::parseError(const char* expected) const
{
    FatalIOError
        << "Parse error, expected a " << expected
        << ", found \n    " << info() << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::token::compound> Foam::token::compound::New
(
    const word& compoundType
)
{
    auto* ctorPtr = emptyConstructorTable(compoundType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "compound",
            compoundType,
            *emptyConstructorTablePtr_
        ) << abort(FatalError);
    }

    return autoPtr<token::compound>(ctorPtr());
}


Foam::autoPtr<Foam::token::compound> Foam::token::compound::New
(
    const word& compoundType,
    Istream& is,
    const bool readContent
)
{
    auto* ctorPtr = emptyConstructorTable(compoundType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            is,
            "compound",
            compoundType,
            *emptyConstructorTablePtr_
        ) << abort(FatalIOError);
    }

    autoPtr<token::compound> aptr(ctorPtr());

    if (readContent)
    {
        aptr->read(is);
    }
    return aptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::token::compound::isCompound(const word& compoundType)
{
    // Could also return the constructor pointer directly
    // and test as bool or use for construction
    //
    // token::compound::emptyConstructorPtr ctorPtr = nullptr;
    // if (emptyConstructorTablePtr_)
    // {
    //     ctorPtr = emptyConstructorTablePtr_->cfind(compoundType);
    //     if (iter.good()) ctorPtr = iter.val();
    // }
    // return ctorPtr;

    return
    (
        emptyConstructorTablePtr_
     && emptyConstructorTablePtr_->contains(compoundType)
    );
}


bool Foam::token::readCompoundToken
(
    const word& compoundType,
    Istream& is,
    const bool readContent
)
{
    // Like compound::New() but more failure tolerant.
    // - a no-op if the compoundType is unknown

    auto* ctorPtr = token::compound::emptyConstructorTable(compoundType);

    if (!ctorPtr)
    {
        return false;
    }

    autoPtr<token::compound> aptr(ctorPtr());

    if (readContent)
    {
        aptr->read(is);
    }
    // Could also set pending(false) for !readContent,
    // but prefer to leave that to the caller.

    (*this) = std::move(aptr);

    return true;
}


Foam::token::compound& Foam::token::transferCompoundToken(const Istream* is)
{
    if (type_ != tokenType::COMPOUND)
    {
        parseError("compound");
    }

    if (data_.compoundPtr->moved())
    {
        if (is)
        {
            FatalIOErrorInFunction(*is)
                << "compound has already been transferred from token\n    "
                << info() << abort(FatalIOError);
        }
        else
        {
            FatalErrorInFunction
                << "compound has already been transferred from token\n    "
                << info() << abort(FatalError);
        }
    }
    // // TDB
    // else if (data_.compoundPtr->pending())
    // {
    //     if (is)
    //     {
    //         FatalIOErrorInFunction(*is)
    //             << "compound is pending (not yet read?)\n    "
    //             << info() << abort(FatalIOError);
    //     }
    //     else
    //     {
    //         FatalErrorInFunction
    //             << "compound is pending (not yet read?)\n    "
    //             << info() << abort(FatalError);
    //     }
    // }
    else
    {
        // TBD: reset pending?
        data_.compoundPtr->moved(true);
    }

    return *data_.compoundPtr;
}


// ************************************************************************* //

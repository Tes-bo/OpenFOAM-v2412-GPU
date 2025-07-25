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

#include "coordSet.H"
#include "coordSetWriter.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

bool Foam::coordSetWriter::supportedType(const word& writeType)
{
    return
    (
        wordConstructorTablePtr_->found(writeType)
     || wordDictConstructorTablePtr_->found(writeType)
    );
}


Foam::autoPtr<Foam::coordSetWriter> Foam::coordSetWriter::New
(
    const word& writeType
)
{
    auto* ctorPtr = wordConstructorTable(writeType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "setWriter",
            writeType,
            *wordConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<coordSetWriter>(ctorPtr());
}


Foam::autoPtr<Foam::coordSetWriter> Foam::coordSetWriter::New
(
    const word& writeType,
    const dictionary& writeOpts
)
{
    // Constructors with dictionary options
    {
        auto* ctorPtr = wordDictConstructorTable(writeType);

        if (ctorPtr)
        {
            return autoPtr<coordSetWriter>(ctorPtr(writeOpts));
        }
    }


    // Constructors without dictionary options
    auto* ctorPtr = wordConstructorTable(writeType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "setWriter",
            writeType,
            *wordConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<coordSetWriter>(ctorPtr());
}


// ************************************************************************* //

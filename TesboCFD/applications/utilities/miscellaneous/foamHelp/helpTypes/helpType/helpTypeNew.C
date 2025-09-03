/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "helpType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::helpType> Foam::helpType::New
(
    const word& helpTypeName
)
{
    auto* ctorPtr = dictionaryConstructorTable(helpTypeName);

    if (!ctorPtr)
    {
        // special treatment for -help
        // exit without stack trace
        if (helpTypeName.starts_with("-help"))
        {
            FatalErrorInFunction
                << "Valid helpType selections:" << nl
                << "    "
                << flatOutput(dictionaryConstructorTablePtr_->sortedToc()) << nl
                << exit(FatalError);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown helpType type '" << helpTypeName << "'" << nl << nl
                << "Valid helpType selections:" << nl
                << "    "
                << flatOutput(dictionaryConstructorTablePtr_->sortedToc()) << nl
                << abort(FatalError);
        }
    }

    Info<< "Selecting helpType '" << helpTypeName << "'" << endl;

    return autoPtr<helpType>(ctorPtr());
}


// ************************************************************************* //

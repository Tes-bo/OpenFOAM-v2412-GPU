/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "edgeMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::edgeMesh> Foam::edgeMesh::New
(
    const fileName& name,
    const word& fileType
)
{
    auto* ctorPtr = fileExtensionConstructorTable(fileType);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Unknown edge format " << fileType
            << " for file " << name << nl << nl
            << "Valid types:" << nl
            << flatOutput(fileExtensionConstructorTablePtr_->sortedToc())
            << exit(FatalError);
    }

    return autoPtr<edgeMesh>(ctorPtr(name));
}


Foam::autoPtr<Foam::edgeMesh> Foam::edgeMesh::New(const fileName& name)
{
    const word ext =
    (
        name.has_ext("gz")
      ? name.stem().ext()
      : name.ext()
    );

    return New(name, ext);
}


// ************************************************************************* //

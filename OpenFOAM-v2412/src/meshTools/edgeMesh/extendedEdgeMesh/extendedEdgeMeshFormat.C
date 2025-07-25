/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "extendedEdgeMeshFormat.H"
#include "IFstream.H"
#include "Time.H"
#include "extendedFeatureEdgeMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::extendedEdgeMeshFormat::extendedEdgeMeshFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::extendedEdgeMeshFormat::read
(
    const fileName& filename
)
{
    clear();

    // Use dummy Time for objectRegistry
    autoPtr<Time> dummyTimePtr(Time::New());

    // Construct IOobject to re-use the headerOk & readHeader
    // (so we can read ascii and binary)
    IOobject io
    (
        filename,
        *dummyTimePtr,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!io.typeHeaderOk<extendedFeatureEdgeMesh>(false))
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    const fileName fName(io.typeFilePath<extendedFeatureEdgeMesh>());

    autoPtr<IFstream> isPtr(new IFstream(fName));
    bool ok = false;
    if (isPtr().good())
    {
        Istream& is = isPtr();
        ok = io.readHeader(is);

        if (ok)
        {
            // Use extendedEdgeMesh IO
            is >> *this;
            ok = is.good();
        }
    }

    return ok;
}


// ************************************************************************* //

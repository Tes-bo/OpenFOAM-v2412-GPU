/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "IOobject.H"
#include "IOstreams.H"
#include "fileOperation.H"  // legacy include
#include "Istream.H"  // legacy include
#include "Pstream.H"  // legacy include

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOobject::typeHeaderOk
(
    const bool checkType,
    const bool search,
    const bool verbose
)
{
    return readAndCheckHeader
    (
        is_globalIOobject<Type>::value,
        Type::typeName,
        checkType,
        search,
        verbose
    );
}


template<class Type>
Foam::fileName Foam::IOobject::typeFilePath(const bool search) const
{
    return
    (
        is_globalIOobject<Type>::value
      ? this->globalFilePath(Type::typeName, search)
      : this->localFilePath(Type::typeName, search)
    );
}


template<class Type>
void Foam::IOobject::warnNoRereading() const
{
    if (readOpt() == IOobjectOption::READ_MODIFIED)
    {
        WarningInFunction
            << Type::typeName << ' ' << name()
            << " constructed with READ_MODIFIED but "
            << Type::typeName << " does not support automatic rereading."
            << endl;
    }
}


// ************************************************************************* //

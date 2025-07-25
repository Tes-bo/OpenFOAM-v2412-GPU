/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "foamGltfObject.H"
#include "foamEndianFwd.H"  // For big/little endian defines

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::glTF::object::object()
:
    base(),
    data_()
{
    #ifdef WM_LITTLE_ENDIAN
    #else
    FatalErrorInFunction
        << "Big-endian buffer support is not available"
        << abort(FatalError);
    #endif
}


Foam::glTF::object::object(const word& name)
:
    base(name),
    data_()
{
    #ifdef WM_LITTLE_ENDIAN
    #else
    FatalErrorInFunction
        << "Big-endian buffer support is not available"
        << abort(FatalError);
    #endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<float>& Foam::glTF::object::data() const noexcept
{
    return data_;
}


// ************************************************************************* //

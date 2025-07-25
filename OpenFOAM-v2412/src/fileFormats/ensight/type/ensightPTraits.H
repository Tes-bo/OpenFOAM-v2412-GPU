/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Class
    Foam::ensightPTraits

Description
    Ensight names and component order for base types.

    For the purpose of traits, integers (label) are treated like
    floating point (scalar). Spherical tensors are mapped as a scalar.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightPTraits_H
#define Foam_ensightPTraits_H

#include "fieldTypes.H"
#include "direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ensightPTraits Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
struct ensightPTraits
{
    //- The type name used in ensight case files
    static const char* const typeName;

    //- Ordering table: return OpenFOAM component given Ensight component
    //  Primarily used for remapping of symmTensor.
    //  OpenFOAM uses (XX, XY, XZ, YY, YZ, ZZ) order,
    //  Ensight uses (XX, YY, ZZ, XY, YZ, XZ) order (like VTK).
    static const direction componentOrder[];
};


// Specializations

template<> const char* const ensightPTraits<label>::typeName;
template<> const direction ensightPTraits<label>::componentOrder[];

template<> const char* const ensightPTraits<float>::typeName;
template<> const direction ensightPTraits<float>::componentOrder[];

template<> const char* const ensightPTraits<double>::typeName;
template<> const direction ensightPTraits<double>::componentOrder[];

template<> const char* const ensightPTraits<vector>::typeName;
template<> const direction ensightPTraits<vector>::componentOrder[];

template<> const char* const ensightPTraits<sphericalTensor>::typeName;
template<> const direction ensightPTraits<sphericalTensor>::componentOrder[];

template<> const char* const ensightPTraits<symmTensor>::typeName;
template<> const direction ensightPTraits<symmTensor>::componentOrder[];

template<> const char* const ensightPTraits<tensor>::typeName;
template<> const direction ensightPTraits<tensor>::componentOrder[];


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

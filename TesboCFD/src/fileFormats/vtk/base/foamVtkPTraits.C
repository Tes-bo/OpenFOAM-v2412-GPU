/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "foamVtkPTraits.H"
#include "foamEndianFwd.H"  // For big/little endian defines

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const
Foam::vtkPTraits<uint8_t>::typeName = "UInt8";

template<>
const char * const
Foam::vtkPTraits<int32_t>::typeName = "Int32";

template<>
const char * const
Foam::vtkPTraits<uint32_t>::typeName = "UInt32";

template<>
const char * const
Foam::vtkPTraits<int64_t>::typeName = "Int64";

template<>
const char * const
Foam::vtkPTraits<uint64_t>::typeName = "UInt64";

template<>
const char * const
Foam::vtkPTraits<float>::typeName = "Float32";

template<>
const char * const
Foam::vtkPTraits<double>::typeName = "Float64";

#ifdef WM_LITTLE_ENDIAN
template<>
const char* const
Foam::vtkPTraits<Foam::endian>::typeName = "LittleEndian";
#else
template<>
const char* const
Foam::vtkPTraits<Foam::endian>::typeName = "BigEndian";
#endif

template<>
const char* const
Foam::vtkPTraits<std::string>::typeName = "String";


// ************************************************************************* //

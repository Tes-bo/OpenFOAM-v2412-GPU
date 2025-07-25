/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

Typedef
    Foam::uLabel

Description
    A uLabel is an uint32_t or uint64_t as specified by the pre-processor macro
    WM_LABEL_SIZE.

    A readULabel function is defined so that uLabel can be constructed from
    Istream.

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_uLabel_H
#define Foam_primitives_uLabel_H

#include "uint.H"
#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UINT_ADD_SIZE(x,s,y) x ## s ## y
#define UINT_ADD_DEF_SIZE(x,s,y) UINT_ADD_SIZE(x,s,y)
#define UINT_SIZE(x,y) UINT_ADD_DEF_SIZE(x,WM_LABEL_SIZE,y)

// Size checks and typedefs (uLabel) in labelFwd.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

constexpr uLabel uLabelMax = UINT_SIZE(UINT, _MAX);

//- Read uLabel from stream.
//  Uses readUint32 or readUint64 according to WM_LABEL_SIZE
inline uLabel readULabel(Istream& is)
{
    return UINT_SIZE(readUint,) (is);
}

//- Parse entire buffer as a uLabel, skipping leading/trailing whitespace.
//  Uses readUint32 or readUint64 according to WM_LABEL_SIZE
//  \return Parsed value or FatalIOError on any problem
inline uLabel readULabel(const char* buf)
{
    return UINT_SIZE(readUint,) (buf);
}

//- Parse entire string as a uLabel, skipping leading/trailing whitespace.
//  Uses readUint32 or readUint64 according to WM_LABEL_SIZE
//  \return Parsed value or FatalIOError on any problem
inline uLabel readULabel(const std::string& str)
{
    return UINT_SIZE(readUint,) (str);
}

//- Parse entire buffer as a uLabel, skipping leading/trailing whitespace.
//  Uses readUint32 or readUint64 according to WM_LABEL_SIZE
//  \return True if successful.
inline bool readULabel(const char* buf, uLabel& val)
{
    return UINT_SIZE(readUint,) (buf, val);
}


//- Parse entire string as a uLabel, skipping leading/trailing whitespace.
//  Uses readUint32 or readUint64 according to WM_LABEL_SIZE
//  \return True if successful.
inline bool readULabel(const std::string& str, uLabel& val)
{
    return UINT_SIZE(readUint,) (str, val);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Raise one uLabel to the power of another
uLabel pow(uLabel a, uLabel b);

//- Evaluate n! : 0 < n <= 12
uLabel factorial(uLabel n);


//- Non-const access to unsigned-type (has no components)
inline uLabel& setComponent(uLabel& val, const direction) noexcept
{
    return val;
}


//- Return unsigned value (has no components)
inline constexpr uLabel component(const uLabel val, const direction) noexcept
{
    return val;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UINT_ADD_SIZE
#undef UINT_ADD_DEF_SIZE
#undef UINT_SIZE

#endif

// ************************************************************************* //

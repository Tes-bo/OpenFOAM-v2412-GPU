/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::label

Description
    A label is an int32_t or int64_t as specified by the pre-processor macro
    WM_LABEL_SIZE.

    A readLabel function is defined so that label can be constructed from
    Istream.

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_label_H
#define Foam_primitives_label_H

#include "int.H"
#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define INT_ADD_SIZE(x,s,y) x ## s ## y
#define INT_ADD_DEF_SIZE(x,s,y) INT_ADD_SIZE(x,s,y)
#define INT_SIZE(x,y) INT_ADD_DEF_SIZE(x,WM_LABEL_SIZE,y)

// Size checks and typedefs (label) in labelFwd.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

constexpr label labelMin = INT_SIZE(INT, _MIN);
constexpr label labelMax = INT_SIZE(INT, _MAX);

//- Parse entire buffer as a label, skipping leading/trailing whitespace.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
//  \return Parsed value or FatalIOError on any problem
inline label readLabel(const char* buf)
{
    return INT_SIZE(readInt,) (buf);
}

//- Parse entire string as a label, skipping leading/trailing whitespace.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
//  \return Parsed value or FatalIOError on any problem
inline label readLabel(const std::string& str)
{
    return INT_SIZE(readInt,) (str);
}

//- Parse entire buffer as a label, skipping leading/trailing whitespace.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
//  \return True if successful.
inline bool readLabel(const char* buf, label& val)
{
    return INT_SIZE(readInt,) (buf, val);
}

//- Parse entire string as a label, skipping leading/trailing whitespace.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
//  \return True if successful.
inline bool readLabel(const std::string& str, label& val)
{
    return INT_SIZE(readInt,) (str, val);
}


//- Read label from stream.
//  Uses readInt32 or readInt64 according to WM_LABEL_SIZE
inline label readLabel(Istream& is)
{
    return INT_SIZE(readInt,) (is);
}

//- Read raw label from binary stream.
//  \note No internal check for binary vs ascii,
//        the caller knows what they are doing
label readRawLabel(Istream& is);

//- Read raw label(s) from binary stream.
//  \note No internal check for binary vs ascii,
//        the caller knows what they are doing
void readRawLabel(Istream& is, label* data, size_t nElem = 1);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Raise one label to the power of another
label pow(label a, label b);

//- Evaluate n! : 0 < n <= 12
label factorial(label n);

//- Non-const access to integer-type (has no components)
inline label& setComponent(label& val, const direction) noexcept
{
    return val;
}

//- Return integer value (has no components)
inline constexpr label component(const label val, const direction) noexcept
{
    return val;
}


// * * * * * * * * * * * * * * * General Functions * * * * * * * * * * * * * //

//- Compare two values for equality
template<class T>
inline bool equal(const T& a, const T& b)
{
    return (a == b);
}

//- Return value clamped between upper and lower limits.
//  Unlike the std::clamp, which selects between references, this version
//  wraps the min/max free functions for component-wise clamping
template<class T>
inline T clamp(const T& val, const T& lower, const T& upper)
{
    return min(max(val, lower), upper);
}


/*---------------------------------------------------------------------------*\
                          Struct labelOp Declaration
\*---------------------------------------------------------------------------*/

//- Conversion/extraction to label operation
//  Specialization of this shall provide a corresponding \c operator().
template<class> struct labelOp;

//- Convert (likely promote) from int32_t to label
template<>
struct labelOp<int32_t>
{
    constexpr label operator()(const int32_t& val) const noexcept
    {
        return val;
    }
};


//- Convert (possibly truncate) from int64_t to label
template<>
struct labelOp<int64_t>
{
    constexpr label operator()(const int64_t& val) const noexcept
    {
        #if WM_LABEL_SIZE == 32
        return label(val);
        #elif WM_LABEL_SIZE == 64
        return val;
        #endif
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Type conversions (narrowing)

//- Type narrowing from int64_t to int32_t
//  Overflow: silently fix, or raise error?
inline int32_t narrowInt32(const int64_t val)
{
    // Single statement - future constexpr?
    return
    (
        (val < INT32_MIN) ? INT32_MIN
      : (val > INT32_MAX) ? INT32_MAX
      : static_cast<int32_t>(val)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelSpecific.H"

#undef INT_ADD_SIZE
#undef INT_ADD_DEF_SIZE
#undef INT_SIZE

#endif

// ************************************************************************* //

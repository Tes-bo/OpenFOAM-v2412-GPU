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

#include "foamVtkLegacyRawFormatter.H"
#include "foamVtkOutputOptions.H"
#include "foamEndian.H"  // For byte swapping
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::legacyRawFormatter::legacyName_ = "BINARY";

const Foam::vtk::outputOptions
Foam::vtk::legacyRawFormatter::opts_(formatType::LEGACY_BINARY);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::vtk::legacyRawFormatter::write
(
    const char* s,
    std::streamsize n
)
{
    os().write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::legacyRawFormatter::legacyRawFormatter
(
    std::ostream& os
)
:
    vtk::formatter(os)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::legacyRawFormatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::legacyRawFormatter::name() const
{
    return legacyName_;
}


const char* Foam::vtk::legacyRawFormatter::encoding() const
{
    return legacyName_;
}


bool Foam::vtk::legacyRawFormatter::writeSize(const uint64_t)
{
    return false;
}


void Foam::vtk::legacyRawFormatter::write
(
    const uint8_t val
)
{
    // Legacy can only handle 32-bit integers.
    // Nonetheless promote to 'label' (32 or 64 bit) and deal with it later
    label copy(val);
    write(copy);
}


void Foam::vtk::legacyRawFormatter::write(const label val)
{
    // std::cerr<<"label is:" << sizeof(val) << '\n';

    // Not entirely correct: the legacy format only supports 32-bit integers.
    // Either limit size for 64-bit label, or simply do not support for 64-bit.
#ifdef WM_LITTLE_ENDIAN
# if WM_LABEL_SIZE == 32
    uint32_t swapped = endian::swap32(val);
    write(reinterpret_cast<const char*>(&swapped), sizeof(uint32_t));
# elif WM_LABEL_SIZE == 64
    uint64_t swapped = endian::swap64(val);
    write(reinterpret_cast<const char*>(&swapped), sizeof(uint64_t));
#endif
#else
    write(reinterpret_cast<const char*>(&val), sizeof(label));
#endif
}


void Foam::vtk::legacyRawFormatter::write(const float val)
{
    // std::cerr<<"float is:" << sizeof(val) << '\n';

#ifdef WM_LITTLE_ENDIAN
    // De-reference in two stages to avoid the warning
    //     dereferencing type-punned pointer will break strict-aliasing rules
    //     [-Wstrict-aliasing]

    const uint32_t* ptr = reinterpret_cast<const uint32_t*>(&val);
    uint32_t swapped = endian::swap32(*ptr);
    write(reinterpret_cast<const char*>(&swapped), sizeof(uint32_t));
#else
    write(reinterpret_cast<const char*>(&val), sizeof(float));
#endif
}


void Foam::vtk::legacyRawFormatter::write(const double val)
{
    // Legacy cannot support Float64 anyhow.
    // std::cerr<<"write double as float:" << val << '\n';

    // Limit range of double to float conversion
    if (val >= std::numeric_limits<float>::max())
    {
        write(std::numeric_limits<float>::max());
    }
    else if (val <= std::numeric_limits<float>::lowest())
    {
        write(std::numeric_limits<float>::lowest());
    }
    else
    {
        float copy(val);
        write(copy);
    }
}


void Foam::vtk::legacyRawFormatter::flush()
{
    os()<< '\n';
}


// ************************************************************************* //

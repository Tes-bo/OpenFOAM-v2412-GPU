/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "foamVtkAppendRawFormatter.H"
#include "foamVtkOutputOptions.H"
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::appendRawFormatter::name_     = "append";
const char* Foam::vtk::appendRawFormatter::encoding_ = "raw";

const Foam::vtk::outputOptions
Foam::vtk::appendRawFormatter::opts_(formatType::APPEND_BINARY);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::vtk::appendRawFormatter::write
(
    const char* s,
    std::streamsize n
)
{
    os().write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::appendRawFormatter::appendRawFormatter(std::ostream& os)
:
    vtk::formatter(os),
    offset_(0)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::appendRawFormatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::appendRawFormatter::name() const
{
    return name_;
}


const char* Foam::vtk::appendRawFormatter::encoding() const
{
    return encoding_;
}


uint64_t Foam::vtk::appendRawFormatter::offset(const uint64_t numbytes)
{
    uint64_t prev = offset_;

    if (formatter::npos != numbytes)
    {
        offset_ += this->encodedLength(sizeof(uint64_t) + numbytes);
    }
    return prev;
}


bool Foam::vtk::appendRawFormatter::writeSize(const uint64_t numbytes)
{
    write(reinterpret_cast<const char*>(&numbytes), sizeof(uint64_t));
    return true;
}


void Foam::vtk::appendRawFormatter::write(const uint8_t val)
{
    write(reinterpret_cast<const char*>(&val), sizeof(uint8_t));
}


void Foam::vtk::appendRawFormatter::write(const label val)
{
    // std::cerr<<"label:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(label));
}


void Foam::vtk::appendRawFormatter::write(const float val)
{
    // std::cerr<<"float:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(float));
}


void Foam::vtk::appendRawFormatter::write(const double val)
{
    // std::cerr<<"double as float=" << val << '\n';

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


void Foam::vtk::appendRawFormatter::flush()
{/*nop*/}


// ************************************************************************* //

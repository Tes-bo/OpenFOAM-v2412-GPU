/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "foamVtkAsciiFormatter.H"
#include "foamVtkOutputOptions.H"
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::asciiFormatter::name_ = "ascii";

const Foam::vtk::outputOptions
Foam::vtk::asciiFormatter::opts_(formatType::INLINE_ASCII);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::vtk::asciiFormatter::next()
{
    if (pos_ >= itemsPerLine_)
    {
        os()<< '\n';
        pos_ = 0;
    }
    else if (pos_)
    {
        os()<< ' ';
    }
    ++pos_;
}


inline void Foam::vtk::asciiFormatter::done()
{
    if (pos_)
    {
        os()<< '\n';
    }
    pos_ = 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::asciiFormatter::asciiFormatter(std::ostream& os)
:
    vtk::formatter(os),
    pos_(0)
{}


Foam::vtk::asciiFormatter::asciiFormatter
(
    std::ostream& os,
    unsigned prec
)
:
    vtk::formatter(os),
    pos_(0)
{
    os.precision(prec);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::asciiFormatter::~asciiFormatter()
{
    done();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::vtk::outputOptions&
Foam::vtk::asciiFormatter::opts() const
{
    return opts_;
}


const char* Foam::vtk::asciiFormatter::name() const
{
    return name_;
}


const char* Foam::vtk::asciiFormatter::encoding() const
{
    return name_;
}


bool Foam::vtk::asciiFormatter::writeSize(const uint64_t)
{
    return false;
}


void Foam::vtk::asciiFormatter::write(const uint8_t val)
{
    next();
    os()<< int(val);
}


void Foam::vtk::asciiFormatter::write(const label val)
{
    next();
    os()<< val;
}


void Foam::vtk::asciiFormatter::write(const float val)
{
    next();
    os()<< val;
}


void Foam::vtk::asciiFormatter::write(const double val)
{
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


void Foam::vtk::asciiFormatter::flush()
{
    done();
}


std::size_t Foam::vtk::asciiFormatter::encodedLength(std::size_t) const
{
    return 0;
}


// ************************************************************************* //

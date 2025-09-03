/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "foamVtkBase64Layer.H"
#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::vtk::foamVtkBase64Layer::encoding_ = "base64";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::vtk::foamVtkBase64Layer::write
(
    const char* s,
    std::streamsize n
)
{
    Foam::base64Layer::write(s, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::foamVtkBase64Layer::foamVtkBase64Layer(std::ostream& os)
:
    vtk::formatter(os),
    Foam::base64Layer(os)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::foamVtkBase64Layer::~foamVtkBase64Layer()
{
    Foam::base64Layer::close();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const char* Foam::vtk::foamVtkBase64Layer::encoding() const
{
    return encoding_;
}


bool Foam::vtk::foamVtkBase64Layer::writeSize(const uint64_t numbytes)
{
    write(reinterpret_cast<const char*>(&numbytes), sizeof(uint64_t));
    return true;
}


void Foam::vtk::foamVtkBase64Layer::write(const uint8_t val)
{
    Foam::base64Layer::add(val);
}


void Foam::vtk::foamVtkBase64Layer::write(const label val)
{
    // std::cerr<<"label:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(label));
}


void Foam::vtk::foamVtkBase64Layer::write(const float val)
{
    // std::cerr<<"float:" << sizeof(val) << "=" << val << '\n';
    write(reinterpret_cast<const char*>(&val), sizeof(float));
}


void Foam::vtk::foamVtkBase64Layer::write(const double val)
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


void Foam::vtk::foamVtkBase64Layer::flush()
{
    Foam::base64Layer::close();
}


std::size_t Foam::vtk::foamVtkBase64Layer::encodedLength(std::size_t n) const
{
    return Foam::base64Layer::encodedLength(n);
}


// ************************************************************************* //

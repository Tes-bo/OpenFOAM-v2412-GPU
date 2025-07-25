/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "base64Layer.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//! \cond fileScope
//- The characters used for base-64
static const unsigned char base64Chars[64] =
{
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
    'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
    'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
    'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
    'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
    'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
    'w', 'x', 'y', 'z', '0', '1', '2', '3',
    '4', '5', '6', '7', '8', '9', '+', '/'
};
//! \endcond


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

std::size_t Foam::base64Layer::encodedLength(std::size_t n)
{
    return 4 * ((n / 3) + (n % 3 ? 1 : 0));
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

inline unsigned char Foam::base64Layer::encode0() const noexcept
{
    // Top 6 bits of char0
    return base64Chars[((group_[0] & 0xFC) >> 2)];
}

inline unsigned char Foam::base64Layer::encode1() const noexcept
{
    // Bottom 2 bits of char0, Top 4 bits of char1
    return base64Chars[((group_[0] & 0x03) << 4) | ((group_[1] & 0xF0) >> 4)];
}

inline unsigned char Foam::base64Layer::encode2() const noexcept
{
    // Bottom 4 bits of char1, Top 2 bits of char2
    return base64Chars[((group_[1] & 0x0F) << 2) | ((group_[2] & 0xC0) >> 6)];
}

inline unsigned char Foam::base64Layer::encode3() const noexcept
{
    // Bottom 6 bits of char2
    return base64Chars[(group_[2] & 0x3F)];
}


void Foam::base64Layer::add(char c)
{
    group_[groupLen_++] = static_cast<unsigned char>(c);
    if (groupLen_ == 3)
    {
        unsigned char out[4];

        out[0] = encode0();
        out[1] = encode1();
        out[2] = encode2();
        out[3] = encode3();
        os_.write(reinterpret_cast<char*>(out), 4);

        groupLen_ = 0;
    }

    dirty_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::base64Layer::base64Layer(std::ostream& os)
:
    os_(os),
    group_(),
    groupLen_(0),
    dirty_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::base64Layer::~base64Layer()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::base64Layer::write(const char* s, std::streamsize n)
{
    for (std::streamsize i=0; i < n; ++i)
    {
        add(s[i]);
    }
}


void Foam::base64Layer::reset()
{
    groupLen_ = 0;
    dirty_ = false;
}


bool Foam::base64Layer::close()
{
    if (!dirty_)
    {
        return false;
    }

    unsigned char out[4];
    if (groupLen_ == 1)
    {
        group_[1] = 0;

        out[0] = encode0();
        out[1] = encode1();
        out[2] = '=';
        out[3] = '=';
        os_.write(reinterpret_cast<char*>(out), 4);
    }
    else if (groupLen_ == 2)
    {
        group_[2] = 0;

        out[0] = encode0();
        out[1] = encode1();
        out[2] = encode2();
        out[3] = '=';
        os_.write(reinterpret_cast<char*>(out), 4);
    }

    // group-length == 0 (no content)
    // group-length == 3 is not possible, already reset in add()

    groupLen_ = 0;
    dirty_ = false;

    return true;
}


// ************************************************************************* //

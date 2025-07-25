/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "FIRECore.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::fileFormats::FIRECore::fileExt3d
>
Foam::fileFormats::FIRECore::file3dExtensions
({
    { fileExt3d::POLY_ASCII, "fpma" },
    { fileExt3d::POLY_BINARY, "fpmb" },
    { fileExt3d::POLY_ASCII_Z, "fpmaz" },
    { fileExt3d::POLY_BINARY_Z, "fpmbz" },
});


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

Foam::label Foam::fileFormats::FIRECore::readPoints
(
    ISstream& is,
    pointField& points
)
{
    const label n = getFireLabel(is);

    if (n > 0)
    {
        points.setSize(n);

        // read the coordinates
        forAll(points, pointI)
        {
            points[pointI] = getFirePoint(is);
        }
    }
    else
    {
        FatalErrorInFunction
            << "no points in file " << is.name()
            << abort(FatalError);
    }

    return n;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileFormats::FIRECore::fireFileName
(
    const fileName& base,
    const enum fileExt3d ext
)
{
    return base + '.' + file3dExtensions[ext];
}


Foam::label Foam::fileFormats::FIRECore::getFireLabel(ISstream& is)
{
    if (is.format() == IOstreamOption::BINARY)
    {
        fireInt_t ivalue;

        is.stdStream().read
        (
            reinterpret_cast<char *>(&ivalue),
            sizeof(ivalue)
        );
        is.syncState();

        return ivalue;
    }
    else
    {
        return readLabel(is);
    }
}


Foam::point Foam::fileFormats::FIRECore::getFirePoint(ISstream& is)
{
    point pt;

    if (is.format() == IOstreamOption::BINARY)
    {
        fireReal_t coord[3];

        is.stdStream().read
        (
            reinterpret_cast<char *>(&coord),
            sizeof(coord)
        );
        is.syncState();

        pt.x() = coord[0];
        pt.y() = coord[1];
        pt.z() = coord[2];
    }
    else
    {
        pt.x() = readScalar(is);
        pt.y() = readScalar(is);
        pt.z() = readScalar(is);
    }

    return pt;
}


std::string Foam::fileFormats::FIRECore::getFireString(ISstream& is)
{
    std::string str;

    if (is.format() == IOstreamOption::BINARY)
    {
        long len;

        is.stdStream().read
        (
            reinterpret_cast<char *>(&len),
            sizeof(len)
        );

        str.resize(len);

        for (std::size_t pos = 0; pos < str.size(); ++pos)
        {
            is.stdStream().read(&(str[pos]), sizeof(char));
        }

        is.syncState();
    }
    else
    {
        const std::string whitespace(" \t\f\v\n\r");

        string s;

        // use a low-level getline, but this means we must handle
        // blank lines manually
        while (s.empty())
        {
            is.getLine(s);
            if (!s.empty())
            {
                // remove prefix whitespace
                size_t pos = s.find_first_not_of(whitespace);

                if (pos != std::string::npos)
                {
                    s.erase(0, pos);

                    // remove suffix whitespace
                    pos = s.find_last_not_of(whitespace);
                    if (pos != std::string::npos)
                    {
                        s.erase(pos + 1);
                    }
                }

                if (pos == std::string::npos)
                {
                    s.clear();
                }
            }
        }

        str.swap(s);
    }

    return str;
}


void Foam::fileFormats::FIRECore::putFireLabel
(
    OSstream& os,
    const label value
)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        fireInt_t ivalue(value);

        os.stdStream().write
        (
            reinterpret_cast<char const *>(&ivalue),
            sizeof(ivalue)
        );

        os.syncState();
    }
    else
    {
        os  << value;
    }
}


void Foam::fileFormats::FIRECore::putFireLabels
(
    OSstream& os,
    const labelUList& lst
)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        fireInt_t ivalue(lst.size());

        os.stdStream().write
        (
            reinterpret_cast<char const *>(&ivalue),
            sizeof(ivalue)
        );

        forAll(lst, i)
        {
            ivalue = lst[i];

            os.stdStream().write
            (
                reinterpret_cast<char const *>(&ivalue),
                sizeof(ivalue)
            );
        }

        os.syncState();
    }
    else
    {
        os  << ' ' << lst.size();
        forAll(lst, i)
        {
            os  << ' ' << lst[i];
        }
        os  << '\n';
    }
}


void Foam::fileFormats::FIRECore::putFireLabels
(
    OSstream& os,
    const label count,
    const label start
)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        fireInt_t ivalue(count);

        os.stdStream().write
        (
            reinterpret_cast<char const *>(&ivalue),
            sizeof(ivalue)
        );

        ivalue = start;
        for (label i=0; i < count; ++i, ++ivalue)
        {
            os.stdStream().write
            (
                reinterpret_cast<char const *>(&ivalue),
                sizeof(ivalue)
            );
        }

        os.syncState();
    }
    else
    {
        os  << ' ' << count;

        label ivalue = start;
        for (label i = 0; i < count; ++i, ++ivalue)
        {
            os  << ' ' << ivalue;
        }
        os  << '\n';
    }
}


void Foam::fileFormats::FIRECore::putFirePoint
(
    OSstream& os,
    const point& value
)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        fireReal_t fvalue[3];
        fvalue[0] = value.x();
        fvalue[1] = value.y();
        fvalue[2] = value.z();

        os.stdStream().write
        (
            reinterpret_cast<char const *>(&fvalue),
            sizeof(fvalue)
        );

        os.syncState();
    }
    else
    {
        os  << ' '
            << value.x() << ' '
            << value.y() << ' '
            << value.z() << '\n';
    }
}


void Foam::fileFormats::FIRECore::putFireString
(
    OSstream& os,
    const std::string& value
)
{
    if (os.format() == IOstreamOption::BINARY)
    {
        long len(value.size());

        os.stdStream().write
        (
            reinterpret_cast<char const *>(&len),
            sizeof(len)
        );

        os.stdStream().write(value.data(), len);
    }
    else
    {
        // output without surrounding quotes
        os.stdStream() << value << '\n';
    }

    os.syncState();
}


// ************************************************************************* //

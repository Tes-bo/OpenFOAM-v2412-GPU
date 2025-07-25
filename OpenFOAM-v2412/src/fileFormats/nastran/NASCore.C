/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "NASCore.H"
#include "IOmanip.H"
#include "Ostream.H"
#include "parsing.H"
#include "defineDebugSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{
    defineDebugSwitchWithName(NASCore, "nastran.debug", 0);
    registerDebugSwitchWithName
    (
        NASCore, Nastran, "nastran.debug"
    );
}
}


const Foam::Enum
<
    Foam::fileFormats::NASCore::fieldFormat
>
Foam::fileFormats::NASCore::fieldFormatNames
({
    { fieldFormat::SHORT, "short" },
    { fieldFormat::LONG,  "long" },
    { fieldFormat::FREE,  "free" },
});


const Foam::Enum
<
    Foam::fileFormats::NASCore::loadFormat
>
Foam::fileFormats::NASCore::loadFormatNames
({
    { loadFormat::PLOAD2, "PLOAD2" },
    { loadFormat::PLOAD4, "PLOAD4" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static inline void putValue(Ostream& os, const Type& value, const int width)
{
    if (width) os << setw(width);
    os  << value;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::scalar Foam::fileFormats::NASCore::readNasScalar(const std::string& str)
{
    const auto signPos = str.find_last_of("+-");

    if
    (
        signPos == std::string::npos
     || signPos == 0
     || str[signPos-1] == 'E' || str[signPos-1] == 'e'
     || isspace(str[signPos-1])
    )
    {
        // A normal number format
        return readScalar(str);
    }


    // Nastran compact number format.
    // Eg, "1234-2" instead of "1234E-2"

    scalar value = 0;
    int exponent = 0; // Any integer

    if
    (
        readScalar(str.substr(0, signPos), value)   // Mantissa
     && readInt(str.substr(signPos), exponent)      // Exponent (with sign)
    )
    {
        // Note: this does not catch underflow/overflow
        // (especially when scalar is a float)
        value *= ::pow(10, exponent);
    }
    else
    {
        FatalIOErrorInFunction("unknown")
            << parsing::errorNames[parsing::errorType::GENERAL] << str
            << exit(FatalIOError);

        value = 0;
    }

    return value;
}


std::string Foam::fileFormats::NASCore::nextNasField
(
    const std::string& str,
    std::string::size_type& pos,
    const std::string::size_type width,
    const bool free_format
)
{
    const auto beg = pos;
    const auto end = str.find(',', pos);

    if (end == std::string::npos)
    {
        if (free_format)
        {
            // Nothing left
            pos = str.size();
            return str.substr(beg);
        }

        // Fixed format - continue after field width
        pos = beg + width;
        return str.substr(beg, width);
    }
    else
    {
        // Free format - continue after comma
        pos = end + 1;
        return str.substr(beg, (end - beg));
    }
}


void Foam::fileFormats::NASCore::setPrecision
(
    Ostream& os,
    const fieldFormat format
)
{
    os.setf(std::ios_base::scientific);

    // Capitalise the E marker
    os.setf(std::ios_base::uppercase);

    const label offset = 7;

    label prec = 16 - offset;
    switch (format)
    {
        case fieldFormat::SHORT :
        {
            prec = 8 - offset;
            break;
        }

        case fieldFormat::LONG :
        case fieldFormat::FREE :
        {
            prec = 16 - offset;
            break;
        }
    }

    os.precision(prec);
}


Foam::Ostream& Foam::fileFormats::NASCore::writeKeyword
(
    Ostream& os,
    const word& keyword,
    const fieldFormat format
)
{
    os.setf(std::ios_base::left);

    switch (format)
    {
        case fieldFormat::SHORT :
        {
            os  << setw(8) << keyword;
            break;
        }
        case fieldFormat::LONG :
        {
            os  << setw(8) << word(keyword + '*');
            break;
        }
        case fieldFormat::FREE :
        {
            os  << keyword;
            break;
        }
    }

    os.unsetf(std::ios_base::left);

    return os;
}


void Foam::fileFormats::NASCore::writeCoord
(
    Ostream& os,
    const point& p,
    const label pointId,  // zero-based
    const fieldFormat format
)
{
    // Field width (SHORT, LONG formats)
    const int width =
    (
        format == fieldFormat::SHORT ? 8
      : format == fieldFormat::LONG ? 16
      : 0
    );

    // Separator char (FREE format)
    const char sep = (format == fieldFormat::FREE ? ',' : '\0');


    // Fixed short/long formats:
    // 1 GRID
    // 2 ID   : point ID - requires starting index of 1
    // 3 CP   : coordinate system ID                (blank)
    // 4 X1   : point x coordinate
    // 5 X2   : point y coordinate
    // 6 X3   : point z coordinate
    // 7 CD   : coordinate system for displacements (blank)
    // 8 PS   : single point constraints            (blank)
    // 9 SEID : super-element ID

    writeKeyword(os, "GRID", format);
    if (sep) os << sep;

    os.setf(std::ios_base::right);

    // Point ID (from 0-based to 1-based)
    putValue(os, (pointId+1), width);
    if (sep) os << sep;

    // Coordinate system ID (blank)
    putValue(os, "", width);
    if (sep) os << sep;

    putValue(os, p.x(), width);
    if (sep) os << sep;

    putValue(os, p.y(), width);
    if (sep) os << sep;

    if (format == fieldFormat::LONG)
    {
        // Continuation
        os.unsetf(std::ios_base::right);
        os << nl;
        writeKeyword(os, "", format);
        os.setf(std::ios_base::right);
    }

    putValue(os, p.z(), width);
    os << nl;

    os.unsetf(std::ios_base::right);
}


Foam::label Foam::fileFormats::NASCore::faceDecomposition
(
    const UList<point>& points,
    const UList<face>& faces,
    labelList& decompOffsets,
    DynamicList<face>& decompFaces
)
{
    // On-demand face decomposition (triangulation)

    decompOffsets.resize(faces.size()+1);
    decompFaces.clear();

    auto offsetIter = decompOffsets.begin();
    *offsetIter = 0; // The first offset is always zero

    for (const face& f : faces)
    {
        const label n = f.size();

        if (n != 3 && n != 4)
        {
            // Decompose non-tri/quad into tris
            f.triangles(points, decompFaces);
        }

        // The end offset, which is the next begin offset
        *(++offsetIter) = decompFaces.size();
    }

    return decompFaces.size();
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "scalarRange.H"
#include "Switch.H"
#include "MinMax.H"
#include "error.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalarRange Foam::scalarRange::always
(
    scalarRange::rangeTypes::ALWAYS,
    -GREAT,
    GREAT
);


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::scalarRange::parse(const std::string& str, scalarRange& range)
{
    range.reset();

    const auto colon = str.find(':');

    if (colon == std::string::npos)
    {
        // No colon

        // Use Switch to accept none/true/false.
        // Others like (f|n|t|y) and (on|off|yes|no) are not really
        // appropriate, but don't worry about that now.

        if (str.size() >= 4)
        {
            Switch sw = Switch::find(str);

            if (sw.good())
            {
                if (sw)
                {
                    range = scalarRange::always;
                }

                return true; // parsed ok
            }
        }

        // "VALUE"
        scalar val;
        if (readScalar(str, val))
        {
            range = scalarRange(val);
        }
    }
    else if (str[colon+1] == ':')
    {
        // A double colon ("::") is a syntax error
        return false;
    }
    else if (colon == 0)
    {
        // ":MAX"
        scalar val;
        if (readScalar(str.substr(1), val))
        {
            range = scalarRange::le(val);
        }
    }
    else if (colon == str.size()-1)
    {
        // "MIN:"
        scalar val;
        if (readScalar(str.substr(0, colon), val))
        {
            range = scalarRange::ge(val);
        }
    }
    else
    {
        // "MIN:MAX"
        scalar minVal, maxVal;
        if
        (
            readScalar(str.substr(0, colon), minVal)
         && readScalar(str.substr(colon+1), maxVal)
        )
        {
            range = scalarRange(minVal, maxVal);
        }
    }

    return range.good();
}


Foam::scalarRange Foam::scalarRange::parse(const std::string& str)
{
    scalarRange range;

    if (!parse(str, range))
    {
        Info<< "Bad scalar-range while parsing: " << str << endl;
    }

    return range;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarRange::scalarRange(const MinMax<label>& range) noexcept
:
    min_(range.min()),
    max_(range.max()),
    type_(max_ < min_ ? rangeTypes::NONE : rangeTypes::GE_LE)
{}


Foam::scalarRange::scalarRange(const MinMax<scalar>& range) noexcept
:
    min_(range.min()),
    max_(range.max()),
    type_(max_ < min_ ? rangeTypes::NONE : rangeTypes::GE_LE)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scalarRange::print(Ostream& os) const
{
    switch (type_)
    {
        case rangeTypes::EQ:
            os << min_;
            break;

        case rangeTypes::GE:
        case rangeTypes::GT:
            os << min_ << ":Inf";
            break;

        case rangeTypes::LE:
        case rangeTypes::LT:
            os << "-Inf:" << max_;
            break;

        case rangeTypes::GE_LE:
            os << min_ << ':' << max_;
            break;

        case rangeTypes::ALWAYS:
            os << "true";
            break;

        default:
            os << "none";
            break;
    }
}


Foam::Ostream& Foam::operator<<(Ostream& os, const scalarRange& range)
{
    range.print(os);
    return os;
}


// ************************************************************************* //

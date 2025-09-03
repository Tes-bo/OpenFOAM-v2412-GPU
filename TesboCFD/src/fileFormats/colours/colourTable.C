/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "colourTable.H"
#include "colourTools.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::colourTable::interpolationType
>
Foam::colourTable::interpolationTypeNames
({
    { interpolationType::RGB, "rgb" },
    { interpolationType::HSV, "hsv" },
    { interpolationType::DIVERGING, "diverging" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::colourTable::colourTable
(
    const List<Tuple2<scalar, vector>>& values,
    const interpolationType interp
)
:
    table_(values),
    interp_(interp)
{}


Foam::colourTable::colourTable
(
    List<Tuple2<scalar, vector>>&& values,
    const interpolationType interp
)
:
    table_(std::move(values)),
    interp_(interp)
{}


Foam::colourTable::colourTable
(
    const dictionary& dict,
    const interpolationType interp
)
:
    table_(),
    interp_(interp)
{
    dict.readEntry("table", table_);
    interpolationTypeNames.readIfPresent("interpolate", dict, interp_);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::colourTable> Foam::colourTable::New(Istream& is)
{
    return autoPtr<colourTable>::New(dictionary(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::colourTable::value(const scalar x) const
{
    if (x <= 0)
    {
        return table_.first().second();
    }

    if (x >= 1)
    {
        return table_.last().second();
    }


    label idx = findLower
    (
        table_, x, 0,
        [](const pair_type& pr, const scalar& val)
        {
            // Test first element
            return (pr.first() <= val);
        }
    );

    if (idx == -1)
    {
        // Use first element only
        return table_.first().second();
    }
    else if (idx == table_.size()-1)
    {
        // Use last element only
        return table_.last().second();
    }

    const scalar t0 = table_[idx].first();
    const scalar t1 = table_[idx+1].first();

    const scalar s = (x - t0)/(t1 - t0);

    const vector& rgb0 = table_[idx].second();
    const vector& rgb1 = table_[idx+1].second();

    if (interp_ == DIVERGING)
    {
        return colourTools::interpolateDiverging(s, rgb0, rgb1);
    }
    else if (interp_ == HSV)
    {
        return colourTools::interpolateHSV(s, rgb0, rgb1);
    }

    return colourTools::interpolateRGB(s, rgb0, rgb1);
}


Foam::List<Foam::Tuple2<Foam::scalar, Foam::vector>>
Foam::colourTable::table(const label nColours) const
{
    List<Tuple2<scalar, vector>> lut(nColours);

    for (label i=0; i < nColours; ++i)
    {
        const scalar x = scalar(i)/scalar(nColours-1);

        lut[i] = pair_type(x, value(x));
    }

    return lut;
}


Foam::Ostream& Foam::colourTable::writeDict(Ostream& os) const
{
    os.beginBlock();
    os.writeEntry("interpolate", interpolationTypeNames[interp_]);
    os.writeEntry("table", table_);
    os.endBlock();

    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const colourTable& tbl)
{
    tbl.writeDict(os);

    return os;
}


// ************************************************************************* //

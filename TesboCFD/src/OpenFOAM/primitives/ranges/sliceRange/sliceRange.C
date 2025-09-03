/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "sliceRange.H"
#include "token.H"
#include "FixedList.H"
#include "List.H"
#include "Istream.H"
#include "Ostream.H"
#include <numeric>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sliceRange::sliceRange(const FixedList<label,3>& coeffs)
:
    start_(coeffs[0]),
    size_(std::max(label(0),coeffs[1])),   // No negative size
    stride_(std::max(label(0),coeffs[2]))  // No negative stride
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::List<Foam::label> Foam::sliceRange::labels() const
{
    List<label> result(size_);

    if (stride_ > 1)
    {
        std::copy(cbegin(), cend(), result.begin());
    }
    else if (stride_ == 1)
    {
        std::iota(result.begin(), result.end(), start_);
    }
    else
    {
        // stride = 0 (identical values!)
        std::fill_n(result.begin(), result.size(), start_);
    }

    return result;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, sliceRange& range)
{
    label beg, len, stride;

    is.readBegin("sliceRange");
    is >> beg >> len >> stride;
    is.readEnd("sliceRange");

    range = sliceRange(beg, len, stride);

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const sliceRange& range)
{
    os  << token::BEGIN_LIST
        << range.start() << token::SPACE
        << range.size() << token::SPACE
        << range.stride() << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

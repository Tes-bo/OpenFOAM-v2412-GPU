/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::objectMap::objectMap(const label index, const labelUList& master)
:
    index_(index),
    objects_(master)
{}


inline Foam::objectMap::objectMap(const label index, labelList&& master)
:
    index_(index),
    objects_(std::move(master))
{}


inline Foam::objectMap::objectMap(Istream& is)
{
    is.readBegin("objectMap");

    is >> index_ >> objects_;

    is.readEnd("objectMap");
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * //

inline bool Foam::operator==(const objectMap& a, const objectMap& b)
{
    return
    (
        (a.index_ == b.index_) && (a.objects_ == b.objects_)
    );
}


inline bool Foam::operator!=(const objectMap& a, const objectMap& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const objectMap& a)
{
    os  << token::BEGIN_LIST
        << a.index_ << token::SPACE
        << a.objects_
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


inline Foam::Istream& Foam::operator>>(Istream& is, objectMap& a)
{
    is.readBegin("objectMap");
    is >> a.index_ >> a.objects_;
    is.readEnd("objectMap");

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //

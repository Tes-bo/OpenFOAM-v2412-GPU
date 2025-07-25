/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "VectorSpace.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::VectorSpace<Form, Cmpt, Ncmpts>::VectorSpace
(
    Istream& is
)
{
    is.readBegin("VectorSpace");

    for (direction i=0; i<Ncmpts; i++)
    {
        is >> v_[i];
    }

    is.readEnd("VectorSpace");

    is.check(FUNCTION_NAME);
}


template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::word Foam::name
(
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    std::ostringstream buf;

    buf << '(' << vs.v_[0];

    for (direction i=1; i<Ncmpts; ++i)
    {
        buf << ',' << vs.v_[i];
    }

    buf << ')';

    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    is.readBegin("VectorSpace");

    for (direction i=0; i<Ncmpts; i++)
    {
        is >> vs.v_[i];
    }

    is.readEnd("VectorSpace");

    is.check(FUNCTION_NAME);

    return is;
}


template<class Form, class Cmpt, Foam::direction Ncmpts>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const VectorSpace<Form, Cmpt, Ncmpts>& vs
)
{
    os << token::BEGIN_LIST << vs.v_[0];

    for (direction i=1; i<Ncmpts; ++i)
    {
        os << token::SPACE << vs.v_[i];
    }

    os << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

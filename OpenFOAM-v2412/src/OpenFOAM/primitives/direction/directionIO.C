/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "direction.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::direction Foam::readDirection(Istream& is)
{
    direction val(0);
    is >> val;

    return val;
}


Foam::Istream& Foam::operator>>(Istream& is, direction& val)
{
    token t(is);

    if (!t.good())
    {
        FatalIOErrorInFunction(is)
            << "Bad token - could not get uint8/direction"
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        val = direction(t.labelToken());
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Wrong token type - expected label (uint8/direction), found "
            << t.info()
            << exit(FatalIOError);
        is.setBad();
        return is;
    }

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const direction val)
{
    os.write(label(val));
    os.check(FUNCTION_NAME);
    return os;
}


std::ostream& Foam::operator<<(std::ostream& os, const direction val)
{
    os << int(val);
    return os;
}


// ************************************************************************* //

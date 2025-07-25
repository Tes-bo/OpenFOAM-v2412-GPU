/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "eddy.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::eddy::eddy(Istream& is)
:
    patchFaceI_(readLabel(is)),
    position0_(is),
    x_(readScalar(is)),
    sigma_(is),
    alpha_(is),
    Rpg_(is),
    c1_(readScalar(is)),
    dir1_(readLabel(is))
{
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

void Foam::eddy::operator=(const eddy& e)
{
    patchFaceI_ = e.patchFaceI_;
    position0_ = e.position0_;
    x_ = e.x_;
    sigma_ = e.sigma_;
    alpha_ = e.alpha_;
    Rpg_ = e.Rpg_;
    c1_ = e.c1_;
    dir1_ = e.dir1_;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, eddy& e)
{
    is.check(FUNCTION_NAME);

    is  >> e.patchFaceI_
        >> e.position0_
        >> e.x_
        >> e.sigma_
        >> e.alpha_
        >> e.Rpg_
        >> e.c1_
        >> e.dir1_;

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const eddy& e)
{
    os.check(FUNCTION_NAME);

    os  << e.patchFaceI_ << token::SPACE
        << e.position0_ << token::SPACE
        << e.x_ << token::SPACE
        << e.sigma_ << token::SPACE
        << e.alpha_ << token::SPACE
        << e.Rpg_ << token::SPACE
        << e.c1_ << token::SPACE
        << e.dir1_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "referredWallFace.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::referredWallFace::referredWallFace
(
    const face& f,
    const pointField& pts,
    label patchi
)
:
    face(f),
    pts_(pts),
    patchi_(patchi)
{
    if (face::size() != pts_.size())
    {
        FatalErrorInFunction
            << "Face and pointField are not the same size." << nl
            << (*this) << nl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::referredWallFace::operator==(const referredWallFace& rhs) const
{
    return
    (
        static_cast<const face&>(rhs) == static_cast<face>(*this)
     && rhs.pts_ == pts_
     && rhs.patchi_ == patchi_
    );
}


bool Foam::referredWallFace::operator!=(const referredWallFace& rhs) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, referredWallFace& rWF)
{
    is  >> static_cast<face&>(rWF) >> rWF.pts_ >> rWF.patchi_;

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const referredWallFace& rWF)
{
    os  << static_cast<const face&>(rWF) << token::SPACE
        << rWF.pts_ << token::SPACE
        << rWF.patchi_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

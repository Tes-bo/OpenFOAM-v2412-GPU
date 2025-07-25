/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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

#include "barycentric2D.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static inline Foam::barycentric2D barycentric2D01Impl
(
    Foam::scalar s,
    Foam::scalar t
)
{
    // Transform the random point in the unit square to a random point in the
    // unit tri by reflecting across the diagonal

    if (s + t > 1)
    {
        s = 1 - s;
        t = 1 - t;
    }

    return Foam::barycentric2D(1 - s - t, s, t);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::barycentric2D Foam::barycentric2D01(Random& rndGen)
{
    const scalar s(rndGen.sample01<scalar>());
    const scalar t(rndGen.sample01<scalar>());

    return barycentric2D01Impl(s, t);
}


// ************************************************************************* //

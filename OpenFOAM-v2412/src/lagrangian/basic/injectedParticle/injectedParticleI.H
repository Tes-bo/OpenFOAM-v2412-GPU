/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::injectedParticle::injectedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    particle(mesh, position, celli),
    position_(position),
    tag_(0),
    soi_(0),
    d_(0),
    U_(Zero)
{}


inline Foam::injectedParticle::injectedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label tag,
    const scalar soi,
    const scalar d,
    const vector& U,
    const bool doLocate
)
:
    particle(mesh, position, -1, -1, -1, doLocate),
    position_(position),
    tag_(tag),
    soi_(soi),
    d_(d),
    U_(U)
{}


// ************************************************************************* //

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

#include "intersection.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::intersection::planarTol_ = 0.2;

const Foam::Enum
<
    Foam::intersection::direction
>
Foam::intersection::directionNames_
({
    { intersection::direction::VECTOR, "vector" },
    { intersection::direction::CONTACT_SPHERE, "contactSphere" },
});


const Foam::Enum
<
    Foam::intersection::algorithm
>
Foam::intersection::algorithmNames_
({
    { intersection::algorithm::FULL_RAY, "fullRay" },
    { intersection::algorithm::HALF_RAY, "halfRay" },
    { intersection::algorithm::VISIBLE, "visible" },
});


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Description
    Correct phi on new faces C-I faces

\*---------------------------------------------------------------------------*/

const surfaceScalarField faceMaskOld
(
    localMin<scalar>(mesh).interpolate(cellMask.oldTime())
);

// Zero Uf on old faceMask (H-I)
rhoUf() *= faceMaskOld;

surfaceVectorField rhoUfint(fvc::interpolate(rho*U));

// Update Uf and phi on new C-I faces
rhoUf() += (1-faceMaskOld)*rhoUfint;

// Update Uf boundary
forAll(rhoUf().boundaryField(), patchI)
{
    rhoUf().boundaryFieldRef()[patchI] = rhoUfint.boundaryField()[patchI];
}

// Calculate absolute flux from the mapped surface velocity
phi = mesh.Sf() & rhoUf();

phi *= faceMask;

// ************************************************************************* //

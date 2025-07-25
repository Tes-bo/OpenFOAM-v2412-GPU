/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "enginePiston.H"
#include "engineTime.H"
#include "polyMesh.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enginePiston::enginePiston
(
    const polyMesh& mesh,
    const word& pistonPatchName,
    const autoPtr<coordinateSystem>& pistonCS,
    const scalar minLayer,
    const scalar maxLayer
)
:
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh.time())),
    patchID_(pistonPatchName, mesh.boundaryMesh()),
    csysPtr_(pistonCS),
    minLayer_(minLayer),
    maxLayer_(maxLayer)
{}


Foam::enginePiston::enginePiston
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh_.time())),
    patchID_(dict.lookup("patch"), mesh.boundaryMesh()),
    csysPtr_
    (
        coordinateSystem::New(mesh_, dict, coordinateSystem::typeName)
    ),
    minLayer_(dict.get<scalar>("minLayer")),
    maxLayer_(dict.get<scalar>("maxLayer"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::enginePiston::writeDict(Ostream& os) const
{
    os  << nl;
    os.beginBlock();
    os.writeEntry("patch", patchID_.name());
    os.writeEntry("minLayer", minLayer_);
    os.writeEntry("maxLayer", maxLayer_);
    os.endBlock();
}


// ************************************************************************* //

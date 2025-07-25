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

#include "slidingInterface.H"
#include "polyTopoChange.H"
#include "polyMesh.H"
#include "polyTopoChanger.H"
#include "polyRemovePoint.H"
#include "polyRemoveFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::clearCouple
(
    polyTopoChange& ref
) const
{
    if (debug)
    {
        Pout<< "void slidingInterface::clearCouple("
            << "polyTopoChange& ref) const for object " << name() << " : "
            << "Clearing old couple points and faces." << endl;
    }

    const polyMesh& mesh = topoChanger().mesh();

    // Remove all points from the point zone
    for
    (
        const label pointi
      : mesh.pointZones()[cutPointZoneID_.index()]
    )
    {
        ref.setAction(polyRemovePoint(pointi));
    }

    // Remove all faces from the face zone
    for
    (
        const label facei
      : mesh.faceZones()[cutFaceZoneID_.index()]
    )
    {
        ref.setAction(polyRemoveFace(facei));
    }

    if (debug)
    {
        Pout<< "void slidingInterface::clearCouple("
            << "polyTopoChange& ref) const for object " << name() << " : "
            << "Finished clearing old couple points and faces." << endl;
    }
}


// ************************************************************************* //

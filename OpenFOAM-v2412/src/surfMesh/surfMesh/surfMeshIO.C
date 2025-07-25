/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "surfMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::setInstance
(
    const fileName& inst,
    IOobjectOption::writeOption wOpt
)
{
    DebugInFunction << "Resetting file instance to " << inst << endl;

    instance() = inst;
    Allocator::setInstance(inst);
    surfZones_.instance()  = inst;

    setWriteOption(wOpt);
}


void Foam::surfMesh::setWriteOption(IOobjectOption::writeOption wOpt)
{
    writeOpt(wOpt);
    Allocator::setWriteOption(wOpt);
    surfZones_.writeOpt(wOpt);
}


Foam::surfMesh::readUpdateState Foam::surfMesh::readUpdate()
{
    DebugInFunction << "Updating mesh based on saved data." << endl;

    // Find point and face instances
    fileName pointsInst(time().findInstance(meshDir(), "points"));
    fileName facesInst(time().findInstance(meshDir(), "faces"));

    DebugInFunction
        << "Points instance: old = " << pointsInstance()
        << " new = " << pointsInst << nl
        << "Faces instance: old = " << facesInstance()
        << " new = " << facesInst << endl;

    if (facesInst != facesInstance())
    {
        // Topological change
        DebugInfo
            << "Topological change" << endl;

        clearOut();

        // Set instance to new instance.
        // Note points instance can differ from faces instance.
        setInstance(facesInst);
        storedIOPoints().instance() = pointsInst;

        storedIOPoints() = pointIOField
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );

        storedFaces() = faceCompactIOList
        (
            IOobject
            (
                "faces",
                facesInst,
                meshSubDir,
                *this,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );

        // Reset the surface zones
        surfZoneIOList newZones
        (
            IOobject
            (
                "surfZones",
                facesInst,
                meshSubDir,
                *this,
                IOobjectOption::LAZY_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );

        // Check that zone types and names are unchanged
        bool zonesChanged = false;

        if (surfZones_.size() != newZones.size())
        {
            zonesChanged = true;
        }
        else
        {
            forAll(surfZones_, zoneI)
            {
                if (surfZones_[zoneI].name() != newZones[zoneI].name())
                {
                    zonesChanged = true;
                    break;
                }
            }
        }

        surfZones_.transfer(newZones);

        if (zonesChanged)
        {
            WarningInFunction
                << "Unexpected consequences.  Proceed with care." << endl;

            return surfMesh::TOPO_PATCH_CHANGE;
        }
        else
        {
            return surfMesh::TOPO_CHANGE;
        }
    }
    else if (pointsInst != pointsInstance())
    {
        // Points moved
        DebugInfo << "Point motion" << endl;

        clearOut();
        storedIOPoints().instance() = pointsInst;

        storedIOPoints() = pointIOField
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );

        return surfMesh::POINTS_MOVED;
    }
    else
    {
        DebugInfo << "No change" << endl;
    }

    return surfMesh::UNCHANGED;
}


bool Foam::surfMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok = Allocator::writeObject(streamOpt, writeOnProc);

    if (ok)
    {
        surfZones_.writeObject(streamOpt, writeOnProc);
    }

    return ok;
}


// ************************************************************************* //

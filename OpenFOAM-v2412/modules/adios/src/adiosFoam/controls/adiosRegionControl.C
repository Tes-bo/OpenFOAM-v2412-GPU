/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "adiosRegionControl.H"
#include "cloud.H"
#include "fvMesh.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

void Foam::adiosFoam::regionControl::read
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    resetAll();
    name_ = mesh.name(); // ensure absolute consistency

    Info<< " Region dict: " << name_ << endl;

    const dictionary* writePtr  = dict.findDict("write");
    const dictionary* ignorePtr = dict.findDict("ignore");

    const dictionary& wrtDict = (writePtr  ? *writePtr  : dictionary::null);
    const dictionary& ignDict = (ignorePtr ? *ignorePtr : dictionary::null);

    explicitWrite_ = wrtDict.getOrDefault("explicit", false);

    requestedFields_.clear();
    requestedClouds_.clear();
    ignoredFields_.clear();
    ignoredClouds_.clear();

    if (explicitWrite_)
    {
        wrtDict.readEntry("fields", requestedFields_);
        wrtDict.readEntry("clouds", requestedClouds_);
    }
    else
    {
        wrtDict.readIfPresent("fields", requestedFields_);
        wrtDict.readIfPresent("clouds", requestedClouds_);
    }

    ignDict.readIfPresent("fields", ignoredFields_);
    ignDict.readIfPresent("clouds", ignoredClouds_);

    // Check if the requested fields are actually accessible
    DynamicList<word> missing(requestedFields_.size());
    for (const wordRe& what : requestedFields_)
    {
        if (what.isLiteral())
        {
            if (mesh.foundObject<regIOobject>(what))
            {
                Info<< " " << what;
            }
            else
            {
                missing.append(what);
            }
        }
    }

    // Also print the cloud names
    for (const auto& name : requestedClouds_)
    {
        Info<< ' ' << name;
    }

    Info<< nl << endl;

    if (missing.size())
    {
        WarningInFunction
            << nl
            << missing.size() << " objects not found in database:" << nl
            << "   " << flatOutput(missing) << endl;
    }
}


void Foam::adiosFoam::regionControl::resetAll()
{
    explicitWrite_ = false;
    requestedFields_.clear();
    requestedClouds_.clear();
    ignoredFields_.clear();
    ignoredClouds_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosFoam::regionControl::regionControl()
:
    regionControl(word::null)
{}


Foam::adiosFoam::regionControl::regionControl(const word& regionName)
:
    name_(regionName),
    explicitWrite_(false),
    requestedFields_(),
    requestedClouds_(),
    ignoredFields_(),
    ignoredClouds_(),
    topoTime_(0),
    pointTime_(0),
    facesName_(),
    pointsName_()
{}


Foam::adiosFoam::regionControl::regionControl
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regionControl(mesh.name())
{
    read(mesh, dict);
    pointTime_ = topoTime_ = mesh.time().timeOutputValue();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::adiosFoam::regionControl::acceptCloudObject
(
    const cloud& c
) const
{
    const word& name = c.name();

    return
    (
        // Not explicitly ignored
        !ignoredClouds_.match(name)
     &&
        (
            // auto-write or explicitly requested
            (autoWrite() && c.writeOpt() == IOobject::AUTO_WRITE)
         || requestedClouds_.match(name)
        )
    );
}


bool Foam::adiosFoam::regionControl::acceptFieldObject
(
    const regIOobject& obj
) const
{
    const word& name = obj.name();

    return
    (
        // NO clouds and not explicitly ignored
        !isA<cloud>(obj) && !ignoredFields_.match(name)
     &&
        (
            // auto-write or explicitly requested
            (autoWrite() && obj.writeOpt() == IOobject::AUTO_WRITE)
         || requestedFields_.match(name)
        )
     && !name.ends_with("_0")  // ignore _0 fields
    );
}


// TODO: handling mesh on input
// void resetPrimitives
// (
// const Xfer<pointField>& points,
// const Xfer<faceList>& faces,
// const Xfer<labelList>& owner,
// const Xfer<labelList>& neighbour,
// const labelList& patchSizes,
// const labelList& patchStarts,
// const bool validBoundary = true
// );
// resetPrimitives(*face,
// Foam::polyMesh::readUpdateState Foam::polyMesh::readUpdate()
// {
//     if (debug)
//     {
//         InfoInFunction << "Updating mesh based on saved data." << endl;
//     }


enum Foam::polyMesh::readUpdateState
Foam::adiosFoam::regionControl::updateTimes(const fvMesh& mesh) const
{
    const double now = mesh.time().timeOutputValue();
    const word pointsInstance = mesh.pointsInstance().name();
    const word facesInstance  = mesh.facesInstance().name();

    polyMesh::readUpdateState state = polyMesh::UNCHANGED;

    if (pointsName_ != pointsInstance)
    {
        state = polyMesh::POINTS_MOVED;

        pointsName_ = pointsInstance;
        pointTime_  = now;
    }

    if (facesName_ != facesInstance)
    {
        state = polyMesh::TOPO_CHANGE;

        facesName_ = facesInstance;
        topoTime_  = now;
    }

    return state;
}


// ************************************************************************* //

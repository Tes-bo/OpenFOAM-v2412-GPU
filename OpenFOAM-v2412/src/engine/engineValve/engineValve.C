/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "engineValve.H"
#include "engineTime.H"
#include "polyMesh.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::engineValve::adjustCrankAngle(const scalar theta) const
{
    if (theta < liftProfileStart_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta < liftProfileStart_)
        {
            adjustedTheta += liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else if (theta > liftProfileEnd_)
    {
        scalar adjustedTheta = theta;

        while (adjustedTheta > liftProfileEnd_)
        {
            adjustedTheta -= liftProfileEnd_ - liftProfileStart_;
        }

        return adjustedTheta;
    }
    else
    {
        return theta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::engineValve::engineValve
(
    const word& name,
    const polyMesh& mesh,
    const autoPtr<coordinateSystem>& valveCS,
    const word& bottomPatchName,
    const word& poppetPatchName,
    const word& stemPatchName,
    const word& curtainInPortPatchName,
    const word& curtainInCylinderPatchName,
    const word& detachInCylinderPatchName,
    const word& detachInPortPatchName,
    const labelList& detachFaces,
    const graph& liftProfile,
    const scalar minLift,
    const scalar minTopLayer,
    const scalar maxTopLayer,
    const scalar minBottomLayer,
    const scalar maxBottomLayer,
    const scalar diameter
)
:
    name_(name),
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh.time())),
    csysPtr_(valveCS.clone()),
    bottomPatch_(bottomPatchName, mesh.boundaryMesh()),
    poppetPatch_(poppetPatchName, mesh.boundaryMesh()),
    stemPatch_(stemPatchName, mesh.boundaryMesh()),
    curtainInPortPatch_(curtainInPortPatchName, mesh.boundaryMesh()),
    curtainInCylinderPatch_(curtainInCylinderPatchName, mesh.boundaryMesh()),
    detachInCylinderPatch_(detachInCylinderPatchName, mesh.boundaryMesh()),
    detachInPortPatch_(detachInPortPatchName, mesh.boundaryMesh()),
    detachFaces_(detachFaces),
    liftProfile_(liftProfile),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(minLift),
    minTopLayer_(minTopLayer),
    maxTopLayer_(maxTopLayer),
    minBottomLayer_(minBottomLayer),
    maxBottomLayer_(maxBottomLayer),
    diameter_(diameter)
{}


Foam::engineValve::engineValve
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    engineDB_(refCast<const engineTime>(mesh_.time())),
    csysPtr_
    (
        coordinateSystem::New(mesh_, dict, coordinateSystem::typeName)
    ),
    bottomPatch_
    (
        dict.get<keyType>("bottomPatch"),
        mesh.boundaryMesh()
    ),
    poppetPatch_
    (
        dict.get<keyType>("poppetPatch"),
        mesh.boundaryMesh()
    ),
    stemPatch_
    (
        dict.get<keyType>("stemPatch"),
        mesh.boundaryMesh()
    ),
    curtainInPortPatch_
    (
        dict.get<keyType>("curtainInPortPatch"),
        mesh.boundaryMesh()
    ),
    curtainInCylinderPatch_
    (
        dict.get<keyType>("curtainInCylinderPatch"),
        mesh.boundaryMesh()
    ),
    detachInCylinderPatch_
    (
        dict.get<keyType>("detachInCylinderPatch"),
        mesh.boundaryMesh()
    ),
    detachInPortPatch_
    (
        dict.get<keyType>("detachInPortPatch"),
        mesh.boundaryMesh()
    ),
    detachFaces_(dict.get<labelList>("detachFaces")),
    liftProfile_("theta", "lift", name_, dict.lookup("liftProfile")),
    liftProfileStart_(min(liftProfile_.x())),
    liftProfileEnd_(max(liftProfile_.x())),
    minLift_(dict.get<scalar>("minLift")),
    minTopLayer_(dict.get<scalar>("minTopLayer")),
    maxTopLayer_(dict.get<scalar>("maxTopLayer")),
    minBottomLayer_(dict.get<scalar>("minBottomLayer")),
    maxBottomLayer_(dict.get<scalar>("maxBottomLayer")),
    diameter_(dict.get<scalar>("diameter"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::engineValve::lift(const scalar theta) const
{
    return interpolateXY
    (
        adjustCrankAngle(theta),
        liftProfile_.x(),
        liftProfile_.y()
    );
}


bool Foam::engineValve::isOpen() const
{
    return lift(engineDB_.theta()) >= minLift_;
}


Foam::scalar Foam::engineValve::curLift() const
{
    return max
    (
        lift(engineDB_.theta()),
        minLift_
    );
}


Foam::scalar Foam::engineValve::curVelocity() const
{
    return
       -(
             curLift()
           - max
             (
                 lift(engineDB_.theta() - engineDB_.deltaTheta()),
                 minLift_
             )
        )/(engineDB_.deltaTValue() + VSMALL);
}


Foam::labelList Foam::engineValve::movingPatchIDs() const
{
    labelList mpIDs(2);
    label nMpIDs = 0;

    if (bottomPatch_.active())
    {
        mpIDs[nMpIDs] = bottomPatch_.index();
        nMpIDs++;
    }

    if (poppetPatch_.active())
    {
        mpIDs[nMpIDs] = poppetPatch_.index();
        nMpIDs++;
    }

    mpIDs.setSize(nMpIDs);

    return mpIDs;
}


void Foam::engineValve::writeDict(Ostream& os) const
{
    os  << nl;
    os.beginBlock(name());

    if (csysPtr_)
    {
        csysPtr_->writeEntry(os);
    }

    os.writeEntry("bottomPatch", bottomPatch_.name());
    os.writeEntry("poppetPatch", poppetPatch_.name());
    os.writeEntry("stemPatch", stemPatch_.name());
    os.writeEntry("curtainInPortPatch", curtainInPortPatch_.name());
    os.writeEntry("curtainInCylinderPatch", curtainInCylinderPatch_.name());
    os.writeEntry("detachInCylinderPatch", detachInCylinderPatch_.name());
    os.writeEntry("detachInPortPatch", detachInPortPatch_.name());
    os.writeEntry("detachFaces", detachFaces_);

    os  << "liftProfile" << nl << token::BEGIN_LIST
        << liftProfile_ << token::END_LIST;
    os.endEntry();

    os.writeEntry("minLift", minLift_);
    os.writeEntry("minTopLayer", minTopLayer_);
    os.writeEntry("maxTopLayer", maxTopLayer_);
    os.writeEntry("minBottomLayer", minBottomLayer_);
    os.writeEntry("maxBottomLayer", maxBottomLayer_);
    os.writeEntry("diameter", diameter_);

    os.endBlock();
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "surfaceDisplacementPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "transformField.H"
#include "fvMesh.H"
#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::surfaceDisplacementPointPatchVectorField::projectMode
>
Foam::surfaceDisplacementPointPatchVectorField::projectModeNames_
({
    { projectMode::NEAREST, "nearest" },
    { projectMode::POINTNORMAL, "pointNormal" },
    { projectMode::FIXEDNORMAL, "fixedNormal" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceDisplacementPointPatchVectorField::calcProjection
(
    vectorField& displacement
) const
{
    const polyMesh& mesh = patch().boundaryMesh().mesh()();
    const labelList& meshPoints = patch().meshPoints();

    //const scalar deltaT = mesh.time().deltaTValue();

    // Construct large enough vector in direction of projectDir so
    // we're guaranteed to hit something.

    //- Per point projection vector:
    const scalar projectLen = mesh.bounds().mag();

    // For case of fixed projection vector:
    vector projectVec(Zero);
    if (projectMode_ == FIXEDNORMAL)
    {
        projectVec = projectLen * normalised(projectDir_);
    }


    // Get fixed points (bit of a hack)
    const pointZone* zonePtr = nullptr;

    if (frozenPointsZone_.size() > 0)
    {
        const pointZoneMesh& pZones = mesh.pointZones();

        zonePtr = &pZones[frozenPointsZone_];

        Pout<< "surfaceDisplacementPointPatchVectorField : Fixing all "
            << zonePtr->size() << " points in pointZone " << zonePtr->name()
            << endl;
    }

    // Get the starting locations from the motionSolver
    const pointField& points0 = mesh.lookupObject<displacementMotionSolver>
    (
        "dynamicMeshDict"
    ).points0();


    pointField start(meshPoints.size());
    forAll(start, i)
    {
        start[i] = points0[meshPoints[i]] + displacement[i];
    }

    label nNotProjected = 0;

    if (projectMode_ == NEAREST)
    {
        List<pointIndexHit> nearest;
        labelList hitSurfaces;
        surfaces().findNearest
        (
            start,
            scalarField(start.size(), sqr(projectLen)),
            hitSurfaces,
            nearest
        );

        forAll(nearest, i)
        {
            const label meshPointi = meshPoints[i];
            const point& pt = mesh.points()[meshPointi];

            if (zonePtr && (zonePtr->whichPoint(meshPointi) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPointi] - pt;
            }
            else if (nearest[i].hit())
            {
                displacement[i] =
                    nearest[i].point()
                  - points0[meshPointi];
            }
            else
            {
                nNotProjected++;

                if (debug)
                {
                    Pout<< "    point:" << meshPointi
                        << " coord:" << pt
                        << "  did not find any surface within " << projectLen
                        << endl;
                }
            }
        }
    }
    else
    {
        // Do tests on all points. Combine later on.

        // 1. Check if already on surface
        List<pointIndexHit> nearest;
        {
            labelList nearestSurface;
            surfaces().findNearest
            (
                start,
                scalarField(start.size(), sqr(SMALL)),
                nearestSurface,
                nearest
            );
        }

        // 2. intersection. (combined later on with information from nearest
        // above)
        vectorField projectVecs(start.size(), projectVec);

        if (projectMode_ == POINTNORMAL)
        {
            projectVecs = projectLen*patch().pointNormals();
        }

        // Knock out any wedge component
        scalarField offset(start.size(), Zero);
        if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
        {
            forAll(offset, i)
            {
                offset[i] = start[i][wedgePlane_];
                start[i][wedgePlane_] = 0;
                projectVecs[i][wedgePlane_] = 0;
            }
        }

        List<pointIndexHit> rightHit;
        {
            labelList rightSurf;
            surfaces().findAnyIntersection
            (
                start,
                start+projectVecs,
                rightSurf,
                rightHit
            );
        }

        List<pointIndexHit> leftHit;
        {
            labelList leftSurf;
            surfaces().findAnyIntersection
            (
                start,
                start-projectVecs,
                leftSurf,
                leftHit
            );
        }

        // 3. Choose either -fixed, nearest, right, left.
        forAll(displacement, i)
        {
            const label meshPointi = meshPoints[i];
            const point& pt = mesh.points()[meshPointi];

            if (zonePtr && (zonePtr->whichPoint(meshPointi) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPointi] - pt;
            }
            else if (nearest[i].hit())
            {
                // Found nearest.
                displacement[i] =
                    nearest[i].point()
                  - points0[meshPointi];
            }
            else
            {
                pointIndexHit interPt;

                if (rightHit[i].hit())
                {
                    if
                    (
                        !leftHit[i].hit()
                    ||
                        (
                            start[i].distSqr(rightHit[i].point())
                          < start[i].distSqr(leftHit[i].point())
                        )
                    )
                    {
                        interPt = rightHit[i];
                    }
                    else
                    {
                        interPt = leftHit[i];
                    }
                }
                else
                {
                    if (leftHit[i].hit())
                    {
                        interPt = leftHit[i];
                    }
                }


                if (interPt.hit())
                {
                    if (wedgePlane_ >= 0 && wedgePlane_ <= vector::nComponents)
                    {
                        interPt.point()[wedgePlane_] += offset[i];
                    }
                    displacement[i] = interPt.point() - points0[meshPointi];
                }
                else
                {
                    nNotProjected++;

                    if (debug)
                    {
                        Pout<< "    point:" << meshPointi
                            << " coord:" << pt
                            << "  did not find any intersection between"
                            << " ray from " << start[i]-projectVecs[i]
                            << " to " << start[i]+projectVecs[i] << endl;
                    }
                }
            }
        }
    }

    reduce(nNotProjected, sumOp<label>());

    if (nNotProjected > 0)
    {
        Info<< "surfaceDisplacement :"
            << " on patch " << patch().name()
            << " did not project " << nNotProjected
            << " out of " << returnReduce(meshPoints.size(), sumOp<label>())
            << " points." << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(p, iF),
    velocity_(Zero),
    projectMode_(NEAREST),
    projectDir_(Zero),
    wedgePlane_(-1)
{}


Foam::surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchVectorField(p, iF, dict),
    velocity_(dict.get<vector>("velocity")),
    surfacesDict_(dict.subDict("geometry")),
    projectMode_(projectModeNames_.get("projectMode", dict)),
    projectDir_
    (
        (projectMode_ == FIXEDNORMAL)
      ? dict.get<vector>("projectDirection")
      : Zero
    ),
    wedgePlane_(dict.getOrDefault("wedgePlane", -1)),
    frozenPointsZone_(dict.getOrDefault("frozenPointsZone", word::null))
{
    if (velocity_.x() < 0 || velocity_.y() < 0 || velocity_.z() < 0)
    {
        FatalErrorInFunction
            << "All components of velocity have to be positive : "
            << velocity_ << nl
            << "Set velocity components to a great value if no clipping"
            << " necessary." << exit(FatalError);
    }
}


Foam::surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const surfaceDisplacementPointPatchVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchVectorField(ppf, p, iF, mapper),
    velocity_(ppf.velocity_),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


Foam::surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const surfaceDisplacementPointPatchVectorField& ppf
)
:
    fixedValuePointPatchVectorField(ppf),
    velocity_(ppf.velocity_),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


Foam::surfaceDisplacementPointPatchVectorField::
surfaceDisplacementPointPatchVectorField
(
    const surfaceDisplacementPointPatchVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchVectorField(ppf, iF),
    velocity_(ppf.velocity_),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::searchableSurfaces&
Foam::surfaceDisplacementPointPatchVectorField::surfaces() const
{
    if (!surfacesPtr_)
    {
        surfacesPtr_.reset
        (
            new searchableSurfaces
            (
                IOobject
                (
                    "abc",                              // dummy name
                    db().time().constant(),             // directory
                    "triSurface",                       // instance
                    db().time(),                        // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfacesDict_,
                true                // allow single-region shortcut
            )
        );
    }

    return *surfacesPtr_;
}


void Foam::surfaceDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const polyMesh& mesh = patch().boundaryMesh().mesh()();

    const vectorField currentDisplacement(this->patchInternalField());

    // Calculate intersections with surface w.r.t points0.
    vectorField displacement(currentDisplacement);
    calcProjection(displacement);


    // offset wrt current displacement
    vectorField offset(displacement-currentDisplacement);

    // Clip offset to maximum displacement possible: velocity*timestep

    const scalar deltaT = mesh.time().deltaTValue();
    const vector clipVelocity = velocity_*deltaT;

    forAll(displacement, i)
    {
        vector& d = offset[i];

        const scalar magD(mag(d));
        if (magD > ROOTVSMALL)
        {
            d /= magD;
            d *= min(magD, mag(clipVelocity));
        }
    }

    this->operator==(currentDisplacement+offset);
    fixedValuePointPatchVectorField::updateCoeffs();
}


void Foam::surfaceDisplacementPointPatchVectorField::write(Ostream& os) const
{
    fixedValuePointPatchField<vector>::write(os);
    os.writeEntry("velocity", velocity_);
    os.writeEntry("geometry", surfacesDict_);
    os.writeEntry("projectMode", projectModeNames_[projectMode_]);
    os.writeEntryIfDifferent<vector>
    (
        "projectDirection",
        Zero,
        projectDir_
    );
    os.writeEntryIfDifferent<label>("wedgePlane", -1, wedgePlane_);
    os.writeEntryIfDifferent<word>
    (
        "frozenPointsZone",
        word::null,
        frozenPointsZone_
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePointPatchTypeField
(
    fixedValuePointPatchVectorField,
    surfaceDisplacementPointPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //

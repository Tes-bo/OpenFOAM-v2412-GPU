/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022,2024 OpenCFD Ltd.
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

#include "surfaceSlipDisplacementPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "transformField.H"
#include "fvMesh.H"
#include "displacementMotionSolver.H"
#include "facePointPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::surfaceSlipDisplacementPointPatchVectorField::projectMode
>
Foam::surfaceSlipDisplacementPointPatchVectorField::projectModeNames_
({
    { projectMode::NEAREST, "nearest" },
    { projectMode::POINTNORMAL, "pointNormal" },
    { projectMode::FIXEDNORMAL, "fixedNormal" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceSlipDisplacementPointPatchVectorField::calcProjection
(
    vectorField& displacement
) const
{
    const polyMesh& mesh = patch().boundaryMesh().mesh()();
    const pointField& localPoints = patch().localPoints();
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

        Pout<< "surfaceSlipDisplacementPointPatchVectorField : Fixing all "
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
            if (zonePtr && (zonePtr->whichPoint(meshPoints[i]) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPoints[i]] - localPoints[i];
            }
            else if (nearest[i].hit())
            {
                displacement[i] = nearest[i].point() - points0[meshPoints[i]];
            }
            else
            {
                nNotProjected++;

                if (debug)
                {
                    Pout<< "    point:" << meshPoints[i]
                        << " coord:" << localPoints[i]
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
        if (wedgePlane_ >= 0 && wedgePlane_ < vector::nComponents)
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
            if (zonePtr && (zonePtr->whichPoint(meshPoints[i]) >= 0))
            {
                // Fixed point. Reset to point0 location.
                displacement[i] = points0[meshPoints[i]] - localPoints[i];
            }
            else if (nearest[i].hit())
            {
                // Found nearest.
                displacement[i] = nearest[i].point() - points0[meshPoints[i]];
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
                    if (wedgePlane_ >= 0 && wedgePlane_ < vector::nComponents)
                    {
                        interPt.point()[wedgePlane_] += offset[i];
                    }
                    displacement[i] = interPt.point() - points0[meshPoints[i]];
                }
                else
                {
                    nNotProjected++;

                    if (debug)
                    {
                        Pout<< "    point:" << meshPoints[i]
                            << " coord:" << localPoints[i]
                            << "  did not find any intersection between"
                            << " ray from " << start[i]-projectVecs[i]
                            << " to " << start[i]+projectVecs[i] << endl;
                    }
                }
            }
        }
    }

    if (scalePtr_)
    {
        const scalarField s
        (
            scalePtr_->value(this->db().time().timeOutputValue())
        );

        displacement *= s;
    }

    reduce(nNotProjected, sumOp<label>());

    if (nNotProjected > 0)
    {
        Info<< "surfaceSlipDisplacement :"
            << " on patch " << patch().name()
            << " did not project " << nNotProjected
            << " out of " << returnReduce(localPoints.size(), sumOp<label>())
            << " points." << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(p, iF),
    projectMode_(NEAREST),
    projectDir_(Zero),
    wedgePlane_(-1)
{}


Foam::surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchVectorField(p, iF, dict),
    surfacesDict_(dict.subDict("geometry")),
    projectMode_(projectModeNames_.get("projectMode", dict)),
    projectDir_
    (
        (projectMode_ == FIXEDNORMAL)
      ? dict.get<vector>("projectDirection")
      : Zero
    ),
    wedgePlane_(dict.getOrDefault("wedgePlane", -1)),
    frozenPointsZone_(dict.getOrDefault("frozenPointsZone", word::null)),
    scalePtr_
    (
        PatchFunction1<scalar>::NewIfPresent
        (
            refCast<const facePointPatch>(p).patch(),
            "scale",
            dict,
            false           // point values
        )
    )
{}


Foam::surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper&
)
:
    pointPatchVectorField(p, iF),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_),
    scalePtr_(ppf.scalePtr_.clone(refCast<const facePointPatch>(p).patch()))
{}


Foam::surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf
)
:
    pointPatchVectorField(ppf),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_),
    scalePtr_
    (
        ppf.scalePtr_.clone
        (
            refCast<const facePointPatch>(ppf.patch()).patch()
        )
    )
{}


Foam::surfaceSlipDisplacementPointPatchVectorField::
surfaceSlipDisplacementPointPatchVectorField
(
    const surfaceSlipDisplacementPointPatchVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(ppf, iF),
    surfacesDict_(ppf.surfacesDict_),
    projectMode_(ppf.projectMode_),
    projectDir_(ppf.projectDir_),
    wedgePlane_(ppf.wedgePlane_),
    frozenPointsZone_(ppf.frozenPointsZone_),
    scalePtr_
    (
        ppf.scalePtr_.clone
        (
            refCast<const facePointPatch>(ppf.patch()).patch()
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::searchableSurfaces&
Foam::surfaceSlipDisplacementPointPatchVectorField::surfaces() const
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
                true    // use single region naming shortcut
            )
        );
    }

    return *surfacesPtr_;
}


void Foam::surfaceSlipDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField displacement(this->patchInternalField());

    // Calculate displacement to project points onto surface
    calcProjection(displacement);

    if (debug)
    {
        Pout<< type() << " :"
            << " on patch " << patch().name()
            << " of field " << this->internalField().name()
            << " projection"
            << " min:" << gMin(displacement)
            << " max:" << gMaxMagSqr(displacement)
            << " average:" << gAverage(displacement)
            << endl;
    }

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->primitiveField());

    //setInInternalField(iF, motionU);
    setInInternalField(iF, displacement);

    pointPatchVectorField::updateCoeffs();
}


void Foam::surfaceSlipDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
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

    if (scalePtr_)
    {
        scalePtr_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePointPatchTypeField
(
    pointPatchVectorField,
    surfaceSlipDisplacementPointPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //

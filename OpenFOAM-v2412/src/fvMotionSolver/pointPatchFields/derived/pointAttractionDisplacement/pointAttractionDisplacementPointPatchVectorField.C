/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "pointAttractionDisplacementPointPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "transformField.H"
#include "displacementMotionSolver.H"
#include "featureEdgeMesh.H"
#include "edgeSlipDisplacementPointPatchVectorField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointAttractionDisplacementPointPatchVectorField::calcProjection
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



    // Get fixed points (bit of a hack)
    const pointZone* zonePtr = nullptr;

    if (frozenPointsZone_.size() > 0)
    {
        const pointZoneMesh& pZones = mesh.pointZones();

        zonePtr = &pZones[frozenPointsZone_];

        Pout<< "pointAttractionDisplacementPointPatchVectorField : Fixing all "
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

    const auto& tree = pointTree();

    label nNotProjected = 0;
    forAll(meshPoints, i)
    {
        const label meshPointi = meshPoints[i];
        const point& pt = mesh.points()[meshPointi];

        if (zonePtr && (zonePtr->whichPoint(meshPointi) >= 0))
        {
            // Fixed point. Reset to point0 location.
            displacement[i] = points0[meshPointi] - pt;
        }
        else
        {
            pointIndexHit nearest = tree.findNearest(start[i], sqr(projectLen));
            if (nearest.hit())
            {
                displacement[i] = nearest.point() - points0[meshPointi];
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

    reduce(nNotProjected, sumOp<label>());

    if (nNotProjected > 0)
    {
        Info<< "pointAttractionDisplacement :"
            << " on patch " << patch().name()
            << " did not project " << nNotProjected
            << " out of " << returnReduce(meshPoints.size(), sumOp<label>())
            << " points." << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointAttractionDisplacementPointPatchVectorField::
pointAttractionDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(p, iF),
    velocity_(Zero)
{}


Foam::pointAttractionDisplacementPointPatchVectorField::
pointAttractionDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    pointPatchVectorField(p, iF, dict),
    velocity_(dict.get<vector>("velocity")),
    featFileName_(dict.get<fileName>("file", keyType::LITERAL)),
    frozenPointsZone_(dict.getOrDefault("frozenPointsZone", word::null))
{
    // Read&store edge mesh on registry
    edgeSlipDisplacementPointPatchVectorField::read
    (
        this->patch().boundaryMesh().mesh().time(),
        dict
    );
}


Foam::pointAttractionDisplacementPointPatchVectorField::
pointAttractionDisplacementPointPatchVectorField
(
    const pointAttractionDisplacementPointPatchVectorField& ppf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper&
)
:
    pointPatchVectorField(p, iF),
    velocity_(ppf.velocity_),
    featFileName_(ppf.featFileName_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


Foam::pointAttractionDisplacementPointPatchVectorField::
pointAttractionDisplacementPointPatchVectorField
(
    const pointAttractionDisplacementPointPatchVectorField& ppf
)
:
    pointPatchVectorField(ppf),
    velocity_(ppf.velocity_),
    featFileName_(ppf.featFileName_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


Foam::pointAttractionDisplacementPointPatchVectorField::
pointAttractionDisplacementPointPatchVectorField
(
    const pointAttractionDisplacementPointPatchVectorField& ppf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    pointPatchVectorField(ppf, iF),
    velocity_(ppf.velocity_),
    featFileName_(ppf.featFileName_),
    frozenPointsZone_(ppf.frozenPointsZone_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataPoint>&
Foam::pointAttractionDisplacementPointPatchVectorField::pointTree() const
{
    if (!pointTreePtr_)
    {
        const Time& tm = this->patch().boundaryMesh().mesh().time();
        const auto& eMesh = tm.lookupObject<edgeMesh>(featFileName_);

        const pointField& points = eMesh.points();

        // Calculate bb of all points
        treeBoundBox bb(points);

        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        bb.inflate(rndGen, 1e-4, ROOTVSMALL);

        pointTreePtr_.reset
        (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint(points),  // All edges

                bb,     // overall search domain
                8,      // maxLevel
                10,     // leafsize
                3.0     // duplicity
            )
        );
    }

    return pointTreePtr_();
}


void Foam::pointAttractionDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const vectorField currentDisplacement(this->patchInternalField());

    // Calculate displacement to project points onto surface
    vectorField displacement(currentDisplacement);
    calcProjection(displacement);


    // offset wrt current displacement
    vectorField offset(displacement-currentDisplacement);

    // Clip offset to maximum displacement possible: velocity*timestep

    const Time& tm = this->patch().boundaryMesh().mesh().time();
    const scalar deltaT = tm.deltaTValue();
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

    // Get internal field to insert values into
    Field<vector>& iF = const_cast<Field<vector>&>(this->primitiveField());

    setInInternalField(iF, (currentDisplacement+offset)());

    pointPatchVectorField::updateCoeffs();
}


void Foam::pointAttractionDisplacementPointPatchVectorField::write
(
    Ostream& os
) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("file", featFileName_);
    os.writeEntryIfDifferent<word>
    (
        "frozenPointsZone",
        word::null,
        frozenPointsZone_
    );
    os.writeEntry("velocity", velocity_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePointPatchTypeField
(
    pointPatchVectorField,
    pointAttractionDisplacementPointPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //

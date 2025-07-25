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

#include "volPointInterpolationAdjoint.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "pointConstraints.H"
#include "surfaceFields.H"
#include "processorPointPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volPointInterpolationAdjoint, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::volPointInterpolationAdjoint::calcBoundaryAddressing()
{
    if (debug)
    {
        Pout<< "volPointInterpolationAdjoint::calcBoundaryAddressing() : "
            << "constructing boundary addressing"
            << endl;
    }

    boundaryPtr_.reset
    (
        new primitivePatch
        (
            SubList<face>
            (
                mesh().faces(),
                mesh().nBoundaryFaces(),
                mesh().nInternalFaces()
            ),
            mesh().points()
        )
    );
    const primitivePatch& boundary = boundaryPtr_();

    boundaryIsPatchFace_.setSize(boundary.size());
    boundaryIsPatchFace_ = false;

    // Store per mesh point whether it is on any 'real' patch. Currently
    // boolList just so we can use syncUntransformedData (does not take
    // bitSet. Tbd)
    boolList isPatchPoint(mesh().nPoints(), false);
    boolList isSymmetryPoint(mesh().nPoints(), false);

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    // Get precalculated volField only so we can use coupled() tests for
    // cyclicAMI
    const surfaceScalarField& magSf = mesh().magSf();

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if
        (
            !isA<emptyPolyPatch>(pp)
         && !magSf.boundaryField()[patchi].coupled()
         && !isA<symmetryPolyPatch>(pp)
         && !isA<symmetryPlanePolyPatch>(pp)
        )
        {
            label bFacei = pp.start()-mesh().nInternalFaces();

            forAll(pp, i)
            {
                boundaryIsPatchFace_[bFacei] = true;

                const face& f = boundary[bFacei++];

                forAll(f, fp)
                {
                    isPatchPoint[f[fp]] = true;
                }
            }
        }
        else if (isA<symmetryPolyPatch>(pp) || isA<symmetryPlanePolyPatch>(pp))
        {
            const labelList& meshPoints = pp.meshPoints();
            for (const label pointI : meshPoints)
            {
                isSymmetryPoint[pointI] = true;
            }
        }
    }

    // Make sure point status is synchronised so even processor that holds
    // no face of a certain patch still can have boundary points marked.
    pointConstraints::syncUntransformedData
    (
        mesh(),
        isPatchPoint,
        orEqOp<bool>()
    );

    // Convert to bitSet
    isPatchPoint_.setSize(mesh().nPoints());
    isPatchPoint_.assign(isPatchPoint);

    isSymmetryPoint_.setSize(mesh().nPoints());
    isSymmetryPoint_.assign(isSymmetryPoint);

    if (debug)
    {
        label nPatchFace = 0;
        forAll(boundaryIsPatchFace_, i)
        {
            if (boundaryIsPatchFace_[i])
            {
                nPatchFace++;
            }
        }
        label nPatchPoint = 0;
        forAll(isPatchPoint_, i)
        {
            if (isPatchPoint_[i])
            {
                nPatchPoint++;
            }
        }
        Pout<< "boundary:" << nl
            << "    faces :" << boundary.size() << nl
            << "    of which on proper patch:" << nPatchFace << nl
            << "    points:" << boundary.nPoints() << nl
            << "    of which on proper patch:" << nPatchPoint << endl;
    }
}


void Foam::volPointInterpolationAdjoint::makeBoundaryWeights
(
    scalarField& sumWeights
)
{
    if (debug)
    {
        Pout<< "volPointInterpolationAdjoint::makeBoundaryWeights() : "
            << "constructing weighting factors for boundary points." << endl;
    }

    const pointField& points = mesh().points();
    const pointField& faceCentres = mesh().faceCentres();

    const primitivePatch& boundary = boundaryPtr_();

    boundaryPointWeights_.clear();
    boundaryPointWeights_.setSize(boundary.meshPoints().size());

    forAll(boundary.meshPoints(), i)
    {
        label pointi = boundary.meshPoints()[i];

        if (isPatchPoint_[pointi])
        {
            const labelList& pFaces = boundary.pointFaces()[i];

            scalarList& pw = boundaryPointWeights_[i];
            pw.setSize(pFaces.size());

            sumWeights[pointi] = 0.0;

            forAll(pFaces, i)
            {
                if (boundaryIsPatchFace_[pFaces[i]])
                {
                    label facei = mesh().nInternalFaces() + pFaces[i];

                    pw[i] = 1.0/mag(points[pointi] - faceCentres[facei]);
                    sumWeights[pointi] += pw[i];
                }
                else
                {
                    pw[i] = 0.0;
                }
            }
        }
    }
}


void Foam::volPointInterpolationAdjoint::makeWeights()
{
    if (debug)
    {
        Pout<< "volPointInterpolationAdjoint::makeWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    const pointMesh& pMesh = pointMesh::New(mesh());

    // Update addressing over all boundary faces
    calcBoundaryAddressing();


    // Running sum of weights
    tmp<pointScalarField> tsumWeights
    (
        new pointScalarField
        (
            IOobject
            (
                "volPointSumWeights",
                mesh().polyMesh::instance(),
                mesh()
            ),
            pMesh,
            dimensionedScalar(dimless, Zero)
        )
    );
    pointScalarField& sumWeights = tsumWeights.ref();


    // Create boundary weights; sumWeights
    makeBoundaryWeights(sumWeights);


    const primitivePatch& boundary = boundaryPtr_();
    const labelList& mp = boundary.meshPoints();


    // Sum collocated contributions
    pointConstraints::syncUntransformedData
    (
        mesh(),
        sumWeights,
        plusEqOp<scalar>()
    );


    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves. Reuse the syncPointData
    // structure.
    pushUntransformedData(sumWeights);

    // Normalise boundary weights
    forAll(mp, i)
    {
        const label pointi = mp[i];

        scalarList& pw = boundaryPointWeights_[i];
        // Note:pw only sized for isPatchPoint
        forAll(pw, i)
        {
            pw[i] /= sumWeights[pointi];
        }
    }

    if (debug)
    {
        Pout<< "volPointInterpolationAdjoint::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::volPointInterpolationAdjoint::volPointInterpolationAdjoint(const fvMesh& vm)
:
    MeshObject_type(vm)
{
    makeWeights();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::volPointInterpolationAdjoint::~volPointInterpolationAdjoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volPointInterpolationAdjoint::updateMesh(const mapPolyMesh&)
{
    makeWeights();
}


bool Foam::volPointInterpolationAdjoint::movePoints()
{
    makeWeights();

    return true;
}


void Foam::volPointInterpolationAdjoint::interpolateSensitivitiesField
(
    const vectorField& pf,
    vectorField& vf,
    const labelHashSet& patchIDs
) const
{
    // Get face data in flat list
    const fvMesh& Mesh = mesh();

    vectorField boundaryVals(Mesh.nBoundaryFaces(), Zero);

    // Do points on 'normal' patches from the surrounding patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const primitivePatch& boundary = boundaryPtr_();
    const labelList& mp = boundary.meshPoints();

    forAll(mp, i)
    {
        label pointi = mp[i];

        if (isPatchPoint_[pointi])
        {
            const labelList& pFaces = boundary.pointFaces()[i];
            const scalarList& pWeights = boundaryPointWeights_[i];
            const vector& val = pf[pointi];
            // In symmetry planes, face-to-point weights should, in general,
            // have half the weight of what they actually do in
            // volPointInterpolation since, in a complete case, a face laying
            // on the opposite side of the symmetry plane would also contribute
            // to a point laying on the symmetry plane.
            // For face-to-point interpolation this is not a problem, but for
            // the adjoint point-to-face interpolation, the correct value of
            // the weight should be taken into consideration
            scalar mod(isSymmetryPoint_[pointi] ? 0.5 : 1);

            forAll(pFaces, j)
            {
                if (boundaryIsPatchFace_[pFaces[j]])
                {
                    boundaryVals[pFaces[j]] += mod*pWeights[j]*val;
                }
            }
        }
    }

    // Transfer values to face-based sensitivity field
    label nPassedFaces(0);
    for (const label patchi : patchIDs)
    {
        const fvPatch& patch = Mesh.boundary()[patchi];
        label bFacei = patch.start() - Mesh.nInternalFaces();
        SubList<vector> patchFaceSens(vf, patch.size(), nPassedFaces);
        patchFaceSens = SubList<vector>(boundaryVals, patch.size(), bFacei);
        nPassedFaces += patch.size();
    }
}


// ************************************************************************* //

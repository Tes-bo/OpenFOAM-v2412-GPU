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

#include "volPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "pointConstraints.H"
#include "surfaceFields.H"
#include "processorPointPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(volPointInterpolation, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::volPointInterpolation::hasSeparated(const pointMesh& pMesh)
{
    const pointBoundaryMesh& pbm = pMesh.boundary();

    bool hasSpecial = false;
    for (const auto& pp : pbm)
    {
        if (isA<coupledFacePointPatch>(pp) && !isType<processorPointPatch>(pp))
        {
            hasSpecial = true;
            break;
        }
    }

    return returnReduceOr(hasSpecial);
}


void Foam::volPointInterpolation::calcBoundaryAddressing()
{
    if (debug)
    {
        Pout<< "volPointInterpolation::calcBoundaryAddressing() : "
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


void Foam::volPointInterpolation::makeInternalWeights(scalarField& sumWeights)
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeInternalWeights() : "
            << "constructing weighting factors for internal and non-coupled"
            << " points." << endl;
    }

    const pointField& points = mesh().points();
    const labelListList& pointCells = mesh().pointCells();
    const vectorField& cellCentres = mesh().cellCentres();

    // Allocate storage for weighting factors
    pointWeights_.clear();
    pointWeights_.setSize(points.size());

    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    forAll(points, pointi)
    {
        if (!isPatchPoint_[pointi])
        {
            const labelList& pcp = pointCells[pointi];

            scalarList& pw = pointWeights_[pointi];
            pw.setSize(pcp.size());

            forAll(pcp, pointCelli)
            {
                pw[pointCelli] =
                    1.0/mag(points[pointi] - cellCentres[pcp[pointCelli]]);

                sumWeights[pointi] += pw[pointCelli];
            }
        }
    }
}


void Foam::volPointInterpolation::makeBoundaryWeights(scalarField& sumWeights)
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeBoundaryWeights() : "
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


void Foam::volPointInterpolation::interpolateOne
(
    const tmp<scalarField>& tnormalisation,
    pointScalarField& pf
) const
{
    if (debug)
    {
        Pout<< "volPointInterpolation::interpolateOne("
            << "pointScalarField&) : "
            << "interpolating oneField"
            <<  " from cells to BOUNDARY points "
            << pf.name() << endl;
    }

    const primitivePatch& boundary = boundaryPtr_();
    const labelList& mp = boundary.meshPoints();
    Field<scalar>& pfi = pf.primitiveFieldRef();


    // 1. Interpolate coupled boundary points from cells
    {
        forAll(mp, i)
        {
            const label pointi = mp[i];
            if (!isPatchPoint_[pointi])
            {
                const scalarList& pw = pointWeights_[pointi];

                scalar& val = pfi[pointi];

                val = Zero;
                forAll(pw, pointCelli)
                {
                    val += pw[pointCelli];
                }
            }
        }
    }


    // 2. Interpolate to the patches preserving fixed value BCs
    {
        forAll(mp, i)
        {
            const label pointi = mp[i];

            if (isPatchPoint_[pointi])
            {
                const labelList& pFaces = boundary.pointFaces()[i];
                const scalarList& pWeights = boundaryPointWeights_[i];

                scalar& val = pfi[pointi];

                val = Zero;
                forAll(pFaces, j)
                {
                    if (boundaryIsPatchFace_[pFaces[j]])
                    {
                        val += pWeights[j];
                    }
                }
            }
        }

        // Sum collocated contributions
        pointConstraints::syncUntransformedData
        (
            mesh(),
            pfi,
            plusEqOp<scalar>()
        );

        // And add separated contributions
        addSeparated(pf);

        // Optionally normalise
        if (tnormalisation)
        {
            const scalarField& normalisation = tnormalisation();
            forAll(mp, i)
            {
                pfi[mp[i]] *= normalisation[i];
            }
        }
    }

    // Apply constraints
    pointConstraints::New(pf.mesh()).constrain(pf, false);
}


void Foam::volPointInterpolation::makeWeights()
{
    if (debug)
    {
        Pout<< "volPointInterpolation::makeWeights() : "
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


    // Create internal weights; add to sumWeights
    makeInternalWeights(sumWeights);


    // Create boundary weights; override sumWeights
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


    // And add separated contributions
    addSeparated(sumWeights);


    // Push master data to slaves. It is possible (not sure how often) for
    // a coupled point to have its master on a different patch so
    // to make sure just push master data to slaves. Reuse the syncPointData
    // structure.
    pushUntransformedData(sumWeights);


    // Normalise internal weights
    forAll(pointWeights_, pointi)
    {
        scalarList& pw = pointWeights_[pointi];
        // Note:pw only sized for !isPatchPoint
        forAll(pw, i)
        {
            pw[i] /= sumWeights[pointi];
        }
    }

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


    // Normalise separated contributions
    if (hasSeparated_)
    {
        if (debug)
        {
            Pout<< "volPointInterpolation::makeWeights() : "
                << "detected separated coupled patches"
                << " - allocating normalisation" << endl;
        }

        // Sum up effect of interpolating one (on boundary points only)
        interpolateOne(tmp<scalarField>(), sumWeights);

        // Store as normalisation factor (on boundary points only)
        normalisationPtr_ = new scalarField(mp.size());
        normalisationPtr_.ref() = scalar(1.0);
        normalisationPtr_.ref() /= scalarField(sumWeights, mp);
    }
    else
    {
        normalisationPtr_.clear();
    }


    if (debug)
    {
        Pout<< "volPointInterpolation::makeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::volPointInterpolation::volPointInterpolation(const fvMesh& vm)
:
    MeshObject_type(vm),
    hasSeparated_(hasSeparated(pointMesh::New(vm)))
{
    makeWeights();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::volPointInterpolation::~volPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::volPointInterpolation::updateMesh(const mapPolyMesh&)
{
    // Recheck whether has coupled patches
    hasSeparated_ = hasSeparated(pointMesh::New(mesh()));

    makeWeights();
}


bool Foam::volPointInterpolation::movePoints()
{
    makeWeights();

    return true;
}


void Foam::volPointInterpolation::interpolateDisplacement
(
    const volVectorField& vf,
    pointVectorField& pf
) const
{
    interpolateInternalField(vf, pf);

    // Interpolate to the patches but no constraints
    interpolateBoundaryField(vf, pf);

    // Apply displacement constraints
    const pointConstraints& pcs = pointConstraints::New(pf.mesh());

    pcs.constrainDisplacement(pf, false);
}


// ************************************************************************* //

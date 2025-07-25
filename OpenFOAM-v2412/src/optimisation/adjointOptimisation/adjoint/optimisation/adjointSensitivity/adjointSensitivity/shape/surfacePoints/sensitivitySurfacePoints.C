/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020, 2022 PCOpt/NTUA
    Copyright (C) 2013-2020, 2022 FOSS GP
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "sensitivitySurfacePoints.H"
#include "deltaBoundary.H"
#include "designVariables.H"
#include "syncTools.H"
#include "symmetryFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivitySurfacePoints, 1);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivitySurfacePoints,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

labelHashSet sensitivitySurfacePoints::populateExtendedIDs() const
{
    // Populate extendedPatchIDs
    label pI(0);
    labelList extendedPatchIDs(mesh_.boundary().size(), -1);
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& pp = mesh_.boundary()[patchI];
        bool isSymmetry
            (isA<symmetryFvPatch>(pp) || isA<symmetryPlaneFvPatch>(pp));
        if (!isA<coupledFvPatch>(pp) && !isA<emptyFvPatch>(pp) && !isSymmetry)
        {
            extendedPatchIDs[pI++] = patchI;
        }
    }
    extendedPatchIDs.setSize(pI);
    return labelHashSet(extendedPatchIDs);
}


void sensitivitySurfacePoints::setSuffixName()
{
    word suffix(adjointMeshMovementSolver_ ? "ESI" : "SI");
    suffix = suffix + word(dict().getOrDefault<word>("suffix", word::null));
    setSuffix(adjointSolver_.solverName() + suffix);
}


void sensitivitySurfacePoints::finalisePointSensitivities()
{
    // List with mesh faces. Global addressing
    const faceList& faces = mesh_.faces();

    // Geometry differentiation engine
    deltaBoundary dBoundary(mesh_);

    for (const label patchI : extendedPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        vectorField nf(patch.nf());

        // Point sens result for patch
        vectorField& pointPatchSens = wallPointSensVecPtr_()[patchI];

        // Face sens for patch
        vectorField facePatchSens = dxdbMult_()[patchI];
        if (dxdbDirectMult_)
        {
            facePatchSens += dxdbDirectMult_()[patchI];
        }
        if (bcDxDbMult_)
        {
            facePatchSens += bcDxDbMult_()[patchI];
        }

        // Correspondance of local point addressing to global point addressing
        const labelList& meshPoints = patch.patch().meshPoints();

        // Each local patch point belongs to these local patch faces
        // (local numbering)
        const labelListList& patchPointFaces = patch.patch().pointFaces();

        // Index of first face in patch
        const label patchStartIndex = patch.start();

        // Loop over patch points.
        // Collect contributions from each boundary face this point belongs to
        forAll(meshPoints, ppI)
        {
            const labelList& pointFaces = patchPointFaces[ppI];
            forAll(pointFaces, pfI)
            {
                label localFaceIndex = pointFaces[pfI];
                label globalFaceIndex = patchStartIndex + localFaceIndex;
                const face& faceI = faces[globalFaceIndex];

                // Point coordinates. All indices in global numbering
                pointField p(faceI.points(mesh_.points()));
                tensorField p_d(faceI.size(), Zero);
                forAll(faceI, facePointI)
                {
                    if (faceI[facePointI] == meshPoints[ppI])
                    {
                        p_d[facePointI] = tensor::I;
                    }
                }
                tensorField deltaNormals =
                    dBoundary.makeFaceCentresAndAreas_d(p, p_d);

                if (isSymmetryPoint_[meshPoints[ppI]])
                {
                    const vector& n = symmPointNormal_[meshPoints[ppI]];
                    deltaNormals =
                      //0.5*(deltaNormals + transform(I - 2.0*sqr(n), deltaNormals));
                        (deltaNormals + transform(I - 2.0*sqr(n), deltaNormals));
                }

                // Element [0] is the variation in the face center
                // (dxFace/dxPoint)
                const tensor& deltaCf = deltaNormals[0];
                pointPatchSens[ppI] += facePatchSens[localFaceIndex] & deltaCf;

                // Term multiplying d(Sf)/d(point displacement) and
                // d(nf)/d(point displacement)
                //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // Element [1] is the variation in the (dimensional) normal
                if (dSfdbMult_)
                {
                    const tensor& deltaSf = deltaNormals[1];
                    pointPatchSens[ppI] +=
                        dSfdbMult_()[patchI][localFaceIndex] & deltaSf;
                }

                // Element [2] is the variation in the unit normal
                if (dnfdbMult_)
                {
                    const tensor& deltaNf = deltaNormals[2];
                    pointPatchSens[ppI] +=
                        dnfdbMult_()[patchI][localFaceIndex] & deltaNf;
                }
            }
        }
    }
}


void sensitivitySurfacePoints::constructGlobalPointNormalsAndAreas
(
    vectorField& pointNormals,
    scalarField& pointMagSf
)
{
    for (const label patchI : extendedPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const scalarField& magSf = patch.magSf();
        vectorField nf(patch.nf());

        // Correspondance of local point addressing to global point addressing
        const labelList& meshPoints = patch.patch().meshPoints();

        // Each local patch point belongs to these local patch faces
        // (local numbering)
        const labelListList& patchPointFaces = patch.patch().pointFaces();

        // Loop over patch points
        forAll(meshPoints, ppI)
        {
            const labelList& pointFaces = patchPointFaces[ppI];
            forAll(pointFaces, pfI)
            {
                const label localFaceIndex = pointFaces[pfI];

                // Accumulate information for point normals
                pointNormals[meshPoints[ppI]] += nf[localFaceIndex];
                pointMagSf[meshPoints[ppI]] += magSf[localFaceIndex];
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pointNormals,
        plusEqOp<vector>(),
        vector::zero
    );
    syncTools::syncPointList
    (
        mesh_,
        pointMagSf,
        plusEqOp<scalar>(),
        scalar(0)
    );

    if (writeGeometricInfo_)
    {
        pointScalarField MagSf
        (
            IOobject
            (
                "pointMagSf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pointMesh::New(mesh_),
            dimensionedScalar(dimless, Zero)
        );
        pointVectorField Nf
        (
            IOobject
            (
                "pointNf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            pointMesh::New(mesh_),
            dimensionedVector(dimless, Zero)
        );
        MagSf.primitiveFieldRef() = pointMagSf;
        Nf.primitiveFieldRef() = pointNormals;
        Nf.primitiveFieldRef().normalise();
        MagSf.write();
        Nf.write();
    }
}


void sensitivitySurfacePoints::computePointDerivativesSize()
{
    // Allocate appropriate space for sensitivities
    label nTotalPoints(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        nTotalPoints += mesh_.boundaryMesh()[patchI].nPoints();
    }

    // Derivatives for all (x,y,z) components of the displacement
    derivatives_ = scalarField(3*nTotalPoints, Zero);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivitySurfacePoints::sensitivitySurfacePoints
(
    const fvMesh& mesh,
    const dictionary& dict,
    adjointSolver& adjointSolver
)
:
    sensitivityShapeESI(mesh, dict, adjointSolver),
    writeGeometricInfo_(false),
    includeSurfaceArea_(false),
    isSymmetryPoint_(mesh.nPoints(), false),
    symmPointNormal_(mesh.nPoints(), Zero),
    extendedPatchIDs_(populateExtendedIDs())
{
    if (debug)
    {
        Info<< "Extended sensitivity patches " << nl;
        for (const label patchI : extendedPatchIDs_)
        {
            Info<< mesh_.boundary()[patchI].name() << endl;
        }
    }
    read();
    setSuffixName();

    // Allocate boundary field pointer
    wallPointSensVecPtr_.reset(createZeroBoundaryPointFieldPtr<vector>(mesh_));
    wallPointSensNormalPtr_.reset
    (
        createZeroBoundaryPointFieldPtr<scalar>(mesh_)
    );
    wallPointSensNormalVecPtr_.reset
    (
        createZeroBoundaryPointFieldPtr<vector>(mesh_)
    );

    computePointDerivativesSize();

    // Populate symmetry patches
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& pp = mesh_.boundary()[patchI];
        bool isSymmetry
            (isA<symmetryFvPatch>(pp) || isA<symmetryPlaneFvPatch>(pp));
        if (isSymmetry)
        {
            const labelList& meshPoints = pp.patch().meshPoints();
            const vectorField& pointNormals = pp.patch().pointNormals();
            forAll(meshPoints, pI)
            {
                const label pointi = meshPoints[pI];
                isSymmetryPoint_[pointi] = true;
                symmPointNormal_[pointi] = pointNormals[pI];
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sensitivitySurfacePoints::read()
{
    writeGeometricInfo_ =
        dict().getOrDefault<bool>("writeGeometricInfo", false);
    // Point sensitivities do not include the surface area by default
    includeSurfaceArea_ =
        dict().getOrDefault<bool>("includeSurfaceArea", false);
}


bool sensitivitySurfacePoints::readDict(const dictionary& dict)
{
    if (sensitivityShapeESI::readDict(dict))
    {
        read();
        return true;
    }

    return false;
}


const Foam::labelHashSet&
Foam::sensitivitySurfacePoints::geometryVariationIntegrationPatches() const
{
    return extendedPatchIDs_;
}


void sensitivitySurfacePoints::assembleSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    // Make sure we have the proper size for the sensitivities
    computePointDerivativesSize();

    // Assemble the multipliers of dxdbFace, as in the ESI approach
    computeDxDbMult();

    // Geometric (or "direct") sensitivities are better computed directly on
    // the points. Compute them and add the ones that depend on dxFace/dxPoint
    finalisePointSensitivities();

    // polyPatch::pointNormals will give the wrong result for points
    // belonging to multiple patches or patch-processorPatch intersections.
    // Keeping a mesh-wide field to allow easy reduction using syncTools.
    // A bit expensive? Better way?
    vectorField pointNormals(mesh_.nPoints(), Zero);
    scalarField pointMagSf(mesh_.nPoints(), Zero);
    constructGlobalPointNormalsAndAreas(pointNormals, pointMagSf);

    // Do parallel communications to avoid wrong values at processor boundaries
    // Global field for accumulation
    vectorField pointSensGlobal(mesh_.nPoints(), Zero);
    for (const label patchI : extendedPatchIDs_)
    {
        const labelList& meshPoints = mesh_.boundaryMesh()[patchI].meshPoints();
        forAll(meshPoints, ppI)
        {
            const label globaPointI = meshPoints[ppI];
            pointSensGlobal[globaPointI] += wallPointSensVecPtr_()[patchI][ppI];
        }
    }

    /*
    // Remove components normal to symmetry planes
    forAll(mesh_.boundary(), patchI)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        if (isA<symmetryFvPatch>(patch) || isA<symmetryPlaneFvPatch>(patch))
        {
            // Deliberately using local point normals instead of the global ones,
            // to get the direction normal to the symmetry plane itself
            const vectorField& pn = patch.patch().pointNormals();
            const labelList& meshPoints = patch.patch().meshPoints();
            forAll(meshPoints, pI)
            {
                const label gpI = meshPoints[pI];
                pointSensGlobal[gpI] -= wallPointSensVecPtr_()[patchI][pI];
                pointSensGlobal[gpI] -=
                    (pointSensGlobal[gpI] & pn[pI])*pn[pI];
                pointSensGlobal[gpI] *= 2;
            }
        }
    }
    */

    // Accumulate dJ/dx_i
    syncTools::syncPointList
    (
        mesh_,
        pointSensGlobal,
        plusEqOp<vector>(),
        vector::zero
    );

    // Transfer back to local fields
    for (const label patchI : extendedPatchIDs_)
    {
        const labelList& meshPoints =
            mesh_.boundaryMesh()[patchI].meshPoints();
        wallPointSensVecPtr_()[patchI].map(pointSensGlobal, meshPoints);
    }

    // Compute normal sens and append to return field
    label nPassedDVs(0);
    const Vector<label>& sd = mesh_.solutionD();
    for (const label patchI : sensitivityPatchIDs_)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        //if (patch.size()>0)
        {
            const labelList& meshPoints = patch.meshPoints();

            // Avoid storing unit point normals in the global list since we
            // might divide multiple times with the number of faces belonging
            // to the point. Instead do the division locally, per patch use
            vectorField patchPointNormals(pointNormals, meshPoints);
            patchPointNormals.normalise();
            if (!includeSurfaceArea_)
            {
                wallPointSensVecPtr_()[patchI] /=
                    scalarField(pointMagSf, meshPoints);
            }
            wallPointSensNormalPtr_()[patchI] =
                wallPointSensVecPtr_()[patchI] & patchPointNormals;
            wallPointSensNormalVecPtr_()[patchI] =
                wallPointSensNormalPtr_()[patchI] *patchPointNormals;

            forAll(patch.localPoints(), pi)
            {
                const label gpi = nPassedDVs + pi;
                const vector& pSens = wallPointSensVecPtr_()[patchI][pi];
                derivatives_[3*gpi    ] = scalar(sd[0] == -1 ? 0 : pSens.x());
                derivatives_[3*gpi + 1] = scalar(sd[1] == -1 ? 0 : pSens.y());
                derivatives_[3*gpi + 2] = scalar(sd[2] == -1 ? 0 : pSens.z());
            }
            nPassedDVs += patch.nPoints();
        }
    }

    // Write derivative fields
    write();

    // Get processed sensitivities from designVariables, if present
    if (designVars)
    {
        adjointSensitivity::assembleSensitivities(designVars);
    }
}


void sensitivitySurfacePoints::write(const word& baseName)
{
    adjointSensitivity::write();
    ShapeSensitivitiesBase::write();

    if (writeGeometricInfo_)
    {
        volVectorField nfOnPatch
        (
            IOobject
            (
                "nfOnPatch",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            Zero
        );

        volVectorField SfOnPatch
        (
            IOobject
            (
                "SfOnPatch",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            Zero
        );

        volVectorField CfOnPatch
        (
            IOobject
            (
                "CfOnPatch",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            Zero
        );
        for (const label patchI : sensitivityPatchIDs_)
        {
            const fvPatch& patch = mesh_.boundary()[patchI];
            nfOnPatch.boundaryFieldRef()[patchI] = patch.nf();
            SfOnPatch.boundaryFieldRef()[patchI] = patch.Sf();
            CfOnPatch.boundaryFieldRef()[patchI] = patch.Cf();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

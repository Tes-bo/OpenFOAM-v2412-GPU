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

#include "sensitivitySurface.H"
#include "volPointInterpolationAdjoint.H"
#include "faMatrices.H"
#include "famSup.H"
#include "famLaplacian.H"
#include "volSurfaceMapping.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sensitivitySurface, 0);
addToRunTimeSelectionTable
(
    adjointSensitivity,
    sensitivitySurface,
    dictionary
);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

label sensitivitySurface::computeFaceDerivativesSize
(
    const bool computeVectorFieldSize
)
{
    label size(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const label patchSize = patch.size();
        size += label(computeVectorFieldSize ? 3*patchSize : patchSize);
    }
    return size;
}


void sensitivitySurface::smoothSensitivities()
{
    // Read in parameters
    const label iters(dict().getOrDefault<label>("iters", 500));
    const scalar tolerance(dict().getOrDefault<scalar>("tolerance", 1.e-06));

    autoPtr<faMesh> aMeshPtr = faMesh::TryNew(mesh_);

    if (aMeshPtr)
    {
        Info<< "Loaded the existing faMesh" << nl;
    }
    else
    {
        // Dictionary used to construct the faMesh
        dictionary faMeshDefinition;

        // Check and read system/faMeshDefinition
        {
            IOobject io
            (
                "faMeshDefinition",
                mesh_.time().caseSystem(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            if (io.typeHeaderOk<IOdictionary>(false))
            {
                Info<< "Using system/faMeshDefinition" << nl;
                faMeshDefinition = IOdictionary(io);
            }
            else if (debug)
            {
                Info<< "No " << io.name() << " in " << io.path() << nl;
            }
        }

        // Check and read system/finite-area/faMeshDefinition
        if (faMeshDefinition.empty())
        {
            IOobject io
            (
                "faMeshDefinition",
                mesh_.time().caseSystem()/faMesh::prefix(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            if (io.typeHeaderOk<IOdictionary>(false))
            {
                Info<< "Using system/finite-area/faMeshDefinition" << nl;
                faMeshDefinition = IOdictionary(io);
            }
            else if (debug)
            {
                Info<< "No " << io.name() << " in " << io.path() << nl;
            }
        }

        // No specified faMeshDefinition?
        // - generate faMesh from all patches on which we compute sensitivities

        if (faMeshDefinition.empty())
        {
            wordList polyMeshPatches(sensitivityPatchIDs_.size());
            label i(0);
            for (const label patchID : sensitivityPatchIDs_)
            {
                polyMeshPatches[i++] = mesh_.boundary()[patchID].name();
            }
            faMeshDefinition.add<wordList>("polyMeshPatches", polyMeshPatches);
            (void)faMeshDefinition.subDictOrAdd("boundary");

            // TBD: Place all edges into the "defaultPatch" ?
            // faMeshDefinition.subDictOrAdd("defaultPatch")
            //     .add("name", "undefined");

            Info<< "Create faMeshDefinition from sensitivity patches"
                << nl << nl;

            faMeshDefinition.writeEntry("faMeshDefinition", Info);
        }

        // Construct faMesh from faMeshDefinition
        aMeshPtr.reset(new faMesh(mesh_, faMeshDefinition));
    }
    faMesh& aMesh = *aMeshPtr;


    // Physical radius of the smoothing, provided either directly or computed
    // based on the average 'length' of boundary faces
    const scalar Rphysical
        (dict().getOrDefault<scalar>("radius", computeRadius(aMesh)));
    DebugInfo
        << "Physical radius of the sensitivity smoothing "
        << Rphysical << nl << endl;

    // Radius used as the diffusivity in the Helmholtz filter, computed as a
    // function of the physical radius
    const dimensionedScalar RpdeSqr
    (
        "RpdeSqr", dimArea, sqr(Rphysical/(2.*::sqrt(3.)))
    );

    dimensionedScalar one(dimless, Foam::one{});

    // Mapping engine
    volSurfaceMapping vsm(aMesh);

    // Source term in faMatrix needs to be an areaField
    areaVectorField sens
    (
        aMesh.newIOobject("sens"),
        aMesh,
        dimensionedVector(dimless, Foam::zero{}),
        faPatchFieldBase::zeroGradientType()
    );

    // Copy sensitivities to area field
    sens.primitiveFieldRef() =
        vsm.mapToSurface<vector>(wallFaceSensVecPtr_());

    // Initialisation of the smoothed sensitivities field based on the original
    // sensitivities
    areaVectorField smoothedSens("smoothedSens", sens);
    for (label iter = 0; iter < iters; ++iter)
    {
        Info<< "Sensitivity smoothing iteration " << iter << endl;

        faVectorMatrix smoothEqn
        (
            fam::Sp(one, smoothedSens)
          - fam::laplacian(RpdeSqr, smoothedSens)
         ==
            sens
        );

        smoothEqn.relax();

        const scalar residual(mag(smoothEqn.solve().initialResidual()));

        DebugInfo
            << "Max smoothSens " << gMax(mag(smoothedSens)()) << endl;

        // Print execution time
        mesh_.time().printExecutionTime(Info);

        // Check convergence
        if (residual < tolerance)
        {
            Info<< "\n***Reached smoothing equation convergence limit, "
                   "iteration " << iter << "***\n\n";
            break;
        }
    }

    // Transfer smooth sensitivity field to wallFaceSensVecPtr_ for defining
    // derivatives_
    vsm.mapToVolume(smoothedSens, wallFaceSensVecPtr_());

    // Write normal, regularised sensitivities to file
    volScalarField volSmoothedSens
    (
        mesh_.newIOobject("smoothedSurfaceSens" + suffix_),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );
    areaVectorField nf(aMesh.faceAreaNormals());
    nf.normalise();
    areaScalarField smoothedSensNormal(smoothedSens & nf);
    vsm.mapToVolume(smoothedSensNormal, volSmoothedSens.boundaryFieldRef());
    volSmoothedSens.write();
}


scalar sensitivitySurface::computeRadius(const faMesh& aMesh)
{
    scalar averageArea(gAverage(aMesh.S().field()));
    const Vector<label>& geometricD = mesh_.geometricD();
    const boundBox& bounds = mesh_.bounds();
    forAll(geometricD, iDir)
    {
        if (geometricD[iDir] == -1)
        {
            averageArea /= bounds.span()[iDir];
        }
    }
    scalar mult = dict().getOrDefault<scalar>("meanRadiusMultiplier", 10);

    return mult*pow(averageArea, scalar(1)/scalar(mesh_.nGeometricD() - 1));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sensitivitySurface::sensitivitySurface
(
    const fvMesh& mesh,
    const dictionary& dict,
    adjointSolver& adjointSolver
)
:
    sensitivitySurfacePoints(mesh, dict, adjointSolver),
    smoothSensitivities_(dict.getOrDefault("smoothSensitivities", false)),
    returnVectorField_
       (dict.getOrDefault<bool>("returnVectorField", true))
    //finalResultIncludesArea_
    //   (dict.getOrDefault<bool>("finalResultIncludesArea", false))
{
    // Allocate boundary field pointers
    wallFaceSensVecPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    wallFaceSensNormalPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    wallFaceSensNormalVecPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));

    derivatives_.setSize(computeFaceDerivativesSize(returnVectorField_), Zero);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void sensitivitySurface::read()
{
    sensitivitySurfacePoints::read();
    smoothSensitivities_ = dict().getOrDefault("smoothSensitivities", false);
    returnVectorField_ =
        dict().getOrDefault<bool>("returnVectorField", true);
    //finalResultIncludesArea_ =
    //   dict().getOrDefault<bool>("finalResultIncludesArea", false);
}


void sensitivitySurface::assembleSensitivities
(
    autoPtr<designVariables>& designVars
)
{
    // Compute point-based sensitivities
    sensitivitySurfacePoints::assembleSensitivities(designVars);

    // Transfer point sensitivities to point field
    vectorField pointSens(mesh_.nPoints(), Zero);
    for (const label patchI : sensitivityPatchIDs_)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];
        const labelList& meshPoints = pp.meshPoints();
        forAll(meshPoints, ppi)
        {
            pointSens[meshPoints[ppi]] = wallPointSensVecPtr_()[patchI][ppi];
        }
    }

    // vectorField face-sensitivities
    vectorField faceVecSens(computeFaceDerivativesSize(false), Zero);

    // Map sensitivities from points to faces
    volPointInterpolationAdjoint interpolation(mesh_);
    interpolation.interpolateSensitivitiesField
        (pointSens, faceVecSens, sensitivityPatchIDs_);

    // Transfer non-regularised sensitivities to wallFaceSens* fields and write
    label nPassedFaces(0);
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        tmp<vectorField> nf = patch.nf();
        wallFaceSensVecPtr_()[patchI] =
            SubField<vector>(faceVecSens, patch.size(), nPassedFaces)
           /patch.magSf();
        wallFaceSensNormalPtr_()[patchI] = wallFaceSensVecPtr_()[patchI] & nf();
        wallFaceSensNormalVecPtr_()[patchI] =
            wallFaceSensNormalPtr_()[patchI]*nf;
        nPassedFaces += patch.size();
    }
    write();

    // Regularise sensitivities if necessary
    if (smoothSensitivities_)
    {
        smoothSensitivities();
    }

    // Make sure we have the correct sensitivities size
    derivatives_.setSize(computeFaceDerivativesSize(returnVectorField_), Zero);
    nPassedFaces = 0;
    for (const label patchI : sensitivityPatchIDs_)
    {
        const fvPatch& patch = mesh_.boundary()[patchI];
        const vectorField nf(patch.nf());
        if (returnVectorField_)
        {
            const Vector<label>& sd = mesh_.solutionD();
            forAll(patch, fI)
            {
                const label gfI = nPassedFaces + fI;
                const vector& fSens = wallFaceSensVecPtr_()[patchI][fI];
                derivatives_[3*gfI    ] = scalar(sd[0] == -1 ? 0 : fSens.x());
                derivatives_[3*gfI + 1] = scalar(sd[1] == -1 ? 0 : fSens.y());
                derivatives_[3*gfI + 2] = scalar(sd[2] == -1 ? 0 : fSens.z());
            }
        }
        else
        {
            forAll(patch, fI)
            {
                derivatives_[nPassedFaces + fI]
                    = wallFaceSensVecPtr_()[patchI][fI] & nf[fI];
            }
        }
        nPassedFaces += patch.size();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

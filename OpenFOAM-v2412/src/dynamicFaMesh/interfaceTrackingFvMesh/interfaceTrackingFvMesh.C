/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "interfaceTrackingFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "motionSolver.H"
#include "volFields.H"
#include "processorFaPatch.H"
#include "wedgeFaPatch.H"
#include "wedgeFaPatchFields.H"
#include "slipFaPatchFields.H"
#include "fixedValueFaPatchFields.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"
#include "wallFvPatch.H"
#include "polyPatchID.H"
#include "fvcMeshPhi.H"
#include "velocityLaplacianFvMotionSolver.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "twoDPointCorrector.H"
#include "gravityMeshObject.H"
#include "turbulentTransportModel.H"
#include "demandDrivenData.H"
#include "unitConversion.H"
#include "foamVtkIndPatchWriter.H"
#include "calculatedFaPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceTrackingFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        interfaceTrackingFvMesh,
        IOobject
    );
    addToRunTimeSelectionTable
    (
        dynamicFvMesh,
        interfaceTrackingFvMesh,
        doInit
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interfaceTrackingFvMesh::initializeData()
{
    // Set free surface patch index
    {
        const word fsPatchName(motion().get<word>("fsPatchName"));

        polyPatchID patch(fsPatchName, this->boundaryMesh());

        if (!patch.active())
        {
            FatalErrorInFunction
                << "Patch name " << fsPatchName << " not found."
                << abort(FatalError);
        }

        fsPatchIndex_ = patch.index();
    }

    // Set point normal correction for finite area mesh
    if (!pointNormalsCorrectionPatches_.empty())
    {
        boolList& correction = aMesh().correctPatchPointNormals();

        for (const word& patchName : pointNormalsCorrectionPatches_)
        {
            label patchID = aMesh().boundary().findPatchID(patchName);

            if (patchID == -1)
            {
                FatalErrorInFunction
                    << "Patch name '" << patchName
                    << "' for point normals correction does not exist"
                    << abort(FatalError);
            }

            correction[patchID] = true;
        }
    }

    // Read motion direction
    if (!normalMotionDir_)
    {
        motionDir_ = normalised(motion().get<vector>("motionDir"));
    }

    // Check if contact angle is defined
    makeContactAngle();

    motion().readIfPresent
    (
        "nonReflectingFreeSurfacePatches",
        nonReflectingFreeSurfacePatches_
    );
}


void Foam::interfaceTrackingFvMesh::makeUs() const
{
    DebugInFunction
        << "making free-surface velocity field" << nl;

    if (UsPtr_)
    {
        FatalErrorInFunction
            << "free-surface velocity field already exists"
            << abort(FatalError);
    }

    wordList patchFieldTypes
    (
        aMesh().boundary().size(),
        faPatchFieldBase::zeroGradientType()
    );

    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == wedgeFaPatch::typeName
        )
        {
            patchFieldTypes[patchI] =
                wedgeFaPatchVectorField::typeName;
        }
        else
        {
            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    mesh().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    patchFieldTypes[patchI] =
                        slipFaPatchVectorField::typeName;
                }
            }
        }
    }

    for (const word& patchName : fixedFreeSurfacePatches_)
    {
        const label fixedPatchID =
            aMesh().boundary().findPatchID(patchName);

        if (fixedPatchID == -1)
        {
            FatalErrorInFunction
                << "Wrong faPatch name '" << patchName
                << "' in the fixedFreeSurfacePatches list"
                << " defined in the dynamicMeshDict dictionary"
                << abort(FatalError);
        }

        label ngbPolyPatchID =
            aMesh().boundary()[fixedPatchID].ngbPolyPatchIndex();

        if (ngbPolyPatchID != -1)
        {
            if
            (
                mesh().boundary()[ngbPolyPatchID].type()
             == wallFvPatch::typeName
            )
            {
                patchFieldTypes[fixedPatchID] =
                    fixedValueFaPatchVectorField::typeName;
            }
        }
    }

    UsPtr_ = new areaVectorField
    (
        IOobject
        (
            "Us",
            aMesh().time().timeName(),
            aMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh(),
        dimensioned<vector>(dimVelocity, Zero),
        patchFieldTypes
    );

    for (const word& patchName : fixedFreeSurfacePatches_)
    {
        const label fixedPatchID = aMesh().boundary().findPatchID(patchName);

        if (fixedPatchID == -1)
        {
            FatalErrorInFunction
                << "Wrong faPatch name '" << patchName
                << "' in the fixedFreeSurfacePatches list"
                << " defined in the dynamicMeshDict dictionary"
                << abort(FatalError);
        }

        label ngbPolyPatchID =
            aMesh().boundary()[fixedPatchID].ngbPolyPatchIndex();

        if (ngbPolyPatchID != -1)
        {
            if
            (
                mesh().boundary()[ngbPolyPatchID].type()
             == wallFvPatch::typeName
            )
            {
                UsPtr_->boundaryFieldRef()[fixedPatchID] == Zero;
            }
        }
    }
}


void Foam::interfaceTrackingFvMesh::makeFsNetPhi() const
{
    DebugInFunction
        << "making free-surface net flux" << nl;

    if (fsNetPhiPtr_)
    {
        FatalErrorInFunction
            << "free-surface net flux already exists"
            << abort(FatalError);
    }

    fsNetPhiPtr_ = new areaScalarField
    (
        IOobject
        (
            "fsNetPhi",
            aMesh().time().timeName(),
            aMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh(),
        dimensionedScalar(dimVelocity*dimArea, Zero)
    );
}


void Foam::interfaceTrackingFvMesh::makeControlPoints()
{
    DebugInFunction
        << "making control points" << nl;

    if (controlPointsPtr_)
    {
        FatalErrorInFunction
            << "control points already exists"
            << abort(FatalError);
    }

    IOobject pointsIO
    (
        "controlPoints",
        mesh().time().timeName(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    );

    if (pointsIO.typeHeaderOk<vectorIOField>())
    {
        Info<< "Reading control points" << endl;
        controlPointsPtr_ = new vectorIOField(pointsIO);
    }
    else
    {
        pointsIO.readOpt(IOobject::NO_READ);

        Info<< "Creating new control points" << endl;
        controlPointsPtr_ = new vectorIOField
        (
            pointsIO,
            aMesh().areaCentres().internalField()
        );

        initializeControlPointsPosition();
    }
}


void Foam::interfaceTrackingFvMesh::makeMotionPointsMask()
{
    DebugInFunction
        << "making motion points mask" << nl;

    if (motionPointsMaskPtr_)
    {
        FatalErrorInFunction
            << "motion points mask already exists"
            << abort(FatalError);
    }

    motionPointsMaskPtr_ = new labelList
    (
        mesh().boundaryMesh()[fsPatchIndex()].nPoints(),
        1
    );

    // Mark free surface boundary points
    // that belong to processor patches
    forAll(aMesh().boundary(), patchI)
    {
        if
        (
            aMesh().boundary()[patchI].type()
         == processorFaPatch::typeName
        )
        {
            const labelList& patchPoints =
                aMesh().boundary()[patchI].pointLabels();

            forAll(patchPoints, pointI)
            {
                motionPointsMask()[patchPoints[pointI]] = -1;
            }
        }
    }

    // Mark fixed free surface boundary points
    for (const word& patchName : fixedFreeSurfacePatches_)
    {
        const label fixedPatchID = aMesh().boundary().findPatchID(patchName);

        if (fixedPatchID == -1)
        {
            FatalErrorInFunction
                << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                << " defined in the dynamicMeshDict dictionary"
                << abort(FatalError);
        }

        const labelList& patchPoints =
            aMesh().boundary()[fixedPatchID].pointLabels();

        forAll(patchPoints, pointI)
        {
            motionPointsMask()[patchPoints[pointI]] = 0;
        }
    }
}


void Foam::interfaceTrackingFvMesh::makeDirections()
{
    DebugInFunction
        << "make displacement directions for points and control points" << nl;

    if (pointsDisplacementDirPtr_ || facesDisplacementDirPtr_)
    {
        FatalErrorInFunction
            << "points, control points displacement directions already exist"
            << abort(FatalError);
    }

    pointsDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[fsPatchIndex()].nPoints(),
            Zero
        );

    facesDisplacementDirPtr_ =
        new vectorField
        (
            mesh().boundaryMesh()[fsPatchIndex()].size(),
            Zero
        );

    if (!normalMotionDir())
    {
        if (mag(motionDir_) < SMALL)
        {
            FatalErrorInFunction
                << "Zero motion direction"
                << abort(FatalError);
        }

        facesDisplacementDir() = motionDir_;
        pointsDisplacementDir() = motionDir_;
    }

    updateDisplacementDirections();
}


void Foam::interfaceTrackingFvMesh::makePhis()
{
    DebugInFunction
        << "making free-surface flux" << nl;

    if (phisPtr_)
    {
        FatalErrorInFunction
            << "free-surface flux already exists"
            << abort(FatalError);
    }

    phisPtr_ = new edgeScalarField
    (
        IOobject
        (
            "phis",
            aMesh().time().timeName(),
            aMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearEdgeInterpolate(Us()) & aMesh().Le()
    );
}


void Foam::interfaceTrackingFvMesh::makeSurfactConc() const
{
    DebugInFunction
        << "making free-surface surfactant concentration field" << nl;

    if (surfactConcPtr_)
    {
        FatalErrorInFunction
            << "free-surface surfactant concentration field already exists"
            << abort(FatalError);
    }

    surfactConcPtr_ = new areaScalarField
    (
        IOobject
        (
            "Cs",
            mesh().time().timeName
            (
                mesh().time().startTime().value()
            ),
            // mesh().time().timeName(),
            aMesh().thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh()
    );
}


void Foam::interfaceTrackingFvMesh::makeBulkSurfactConc() const
{
    DebugInFunction
        << "making volume surfactant concentration field" << nl;

    if (bulkSurfactConcPtr_)
    {
        FatalErrorInFunction
            << "volume surfactant concentration field already exists"
            << abort(FatalError);
    }

    bulkSurfactConcPtr_ = new volScalarField
    (
        IOobject
        (
            "C",
            mesh().time().timeName
            (
                mesh().time().startTime().value()
            ),
            // mesh().time().timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    );
    volScalarField& bulkSurfactConc = *bulkSurfactConcPtr_;

    if (mesh().time().timeIndex()-1 == 0)
    {
        // Initialize uniform volume surfactant concentration
        bulkSurfactConc = surfactant().bulkConc();
        bulkSurfactConc.correctBoundaryConditions();
    }
}


void Foam::interfaceTrackingFvMesh::makeSurfaceTension() const
{
    DebugInFunction
        << "making surface tension field" << nl;

    if (surfaceTensionPtr_)
    {
        FatalErrorInFunction
            << "surface tension field already exists"
            << abort(FatalError);
    }

    surfaceTensionPtr_ = new areaScalarField
    (
        IOobject
        (
            "surfaceTension",
            aMesh().time().timeName(),
            aMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sigma() + surfactant().dSigma(surfactantConcentration())/rho_
    );
}


void Foam::interfaceTrackingFvMesh::makeSurfactant() const
{
    DebugInFunction
        << "making surfactant properties" << nl;

    if (surfactantPtr_)
    {
        FatalErrorInFunction
            << "surfactant properties already exists"
            << abort(FatalError);
    }

    const dictionary& surfactProp =
        motion().subDict("surfactantProperties");

    surfactantPtr_ = new surfactantProperties(surfactProp);
}


void Foam::interfaceTrackingFvMesh::makeContactAngle()
{
    DebugInFunction
        << "making contact angle field" << nl;

    if (contactAnglePtr_)
    {
        FatalErrorInFunction
            << "contact angle already exists"
            << abort(FatalError);
    }

    // Check if contactAngle is defined
    IOobject contactAngleHeader
    (
        "contactAngle",
        aMesh().time().timeName(),
        aMesh().thisDb(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (contactAngleHeader.typeHeaderOk<areaScalarField>())
    {
        Info<< "Reading contact angle field" << endl;

        contactAnglePtr_ = new areaScalarField(contactAngleHeader, aMesh());
    }
}


void Foam::interfaceTrackingFvMesh::updateDisplacementDirections()
{
    if (normalMotionDir())
    {
        // Update point displacement direction
        pointsDisplacementDir() = aMesh().pointAreaNormals();

        // Correct point displacement direction at contact line
        forAll(aMesh().boundary(), patchI)
        {
            if (contactAnglePtr_)
            {
                label ngbPolyPatchID =
                    aMesh().boundary()[patchI].ngbPolyPatchIndex();

                if (ngbPolyPatchID != -1)
                {
                    if
                    (
                        mesh().boundary()[ngbPolyPatchID].type()
                     == wallFvPatch::typeName
                    )
                    {
                        labelList patchPoints =
                            aMesh().boundary()[patchI].pointLabels();

                        vectorField N
                        (
                            aMesh().boundary()[patchI]
                           .ngbPolyPatchPointNormals()
                        );

                        forAll(patchPoints, pointI)
                        {
                            pointsDisplacementDir()[patchPoints[pointI]] -=
                                N[pointI]
                               *(
                                   N[pointI]
                                 & pointsDisplacementDir()[patchPoints[pointI]]
                                );

                            pointsDisplacementDir()[patchPoints[pointI]] /=
                                mag
                                (
                                    pointsDisplacementDir()
                                    [
                                        patchPoints[pointI]
                                    ]
                                ) + SMALL;
                        }
                    }
                }
            }
        }

        // Update face displacement direction
        facesDisplacementDir() =
            aMesh().faceAreaNormals().internalField();

        // Correction of control points position
        const vectorField& Cf = aMesh().areaCentres().internalField();

        controlPoints() =
            facesDisplacementDir()
           *(facesDisplacementDir()&(controlPoints() - Cf))
          + Cf;
    }
}


void Foam::interfaceTrackingFvMesh::initializeControlPointsPosition()
{
    {
        const faceList& faces = aMesh().faces();
        const pointField& points = aMesh().points();

        pointField displacement(pointDisplacement());
        scalarField sweptVolCorr(faces.size(), Zero);
        correctPointDisplacement(sweptVolCorr, displacement);

        pointField newPoints(points + displacement);

        scalarField sweptVol(faces.size(), Zero);

        forAll(faces, faceI)
        {
            sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
        }

        vectorField faceArea(faces.size(), Zero);

        forAll(faceArea, faceI)
        {
            faceArea[faceI] = faces[faceI].unitNormal(newPoints);
        }

        scalarField deltaH = scalarField(aMesh().nFaces(), Zero);

        forAll(deltaH, faceI)
        {
            deltaH[faceI] = sweptVol[faceI]/
                ((faceArea[faceI] & facesDisplacementDir()[faceI]) + SMALL);

            if (mag(faceArea[faceI] & facesDisplacementDir()[faceI]) < SMALL)
            {
                // Info<< (faceArea[faceI] & facesDisplacementDir()[faceI])
                //     << ", " << faceArea[faceI]
                //     << ", " << facesDisplacementDir()[faceI] << endl;

                FatalError
                    << "Something is wrong with specified motion direction"
                    << abort(FatalError);
            }
        }

        for (const word& patchName : fixedFreeSurfacePatches_)
        {
            label fixedPatchID = aMesh().boundary().findPatchID(patchName);

            if (fixedPatchID == -1)
            {
                FatalError
                    << "Wrong faPatch name in the fixedFreeSurfacePatches list"
                    << " defined in the freeSurfaceProperties dictionary"
                    << abort(FatalError);
            }

            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

            forAll(eFaces, edgeI)
            {
                deltaH[eFaces[edgeI]] *= 2.0;
            }
        }

        controlPoints() += facesDisplacementDir()*deltaH;
    }
}


void Foam::interfaceTrackingFvMesh::smoothFreeSurfaceMesh()
{
    Info<< "Smoothing free surface mesh" << endl;

    controlPoints() = aMesh().areaCentres().internalField();

    pointField displacement(pointDisplacement());

    const faceList& faces = aMesh().faces();
    const pointField& points = aMesh().points();

    pointField newPoints(points + displacement);

    scalarField sweptVol(faces.size(), Zero);
    forAll(faces, faceI)
    {
        sweptVol[faceI] = -faces[faceI].sweptVol(points, newPoints);
    }

    vectorField faceArea(faces.size(), Zero);
    forAll(faceArea, faceI)
    {
        faceArea[faceI] = faces[faceI].unitNormal(newPoints);
    }

    scalarField deltaHf
    (
        sweptVol/(faceArea & facesDisplacementDir())
    );

    for (const word& patchName : fixedFreeSurfacePatches_)
    {
        label fixedPatchID = aMesh().boundary().findPatchID(patchName);

        if (fixedPatchID == -1)
        {
            FatalError
                << "Wrong faPatch name fixedFreeSurfacePatches list"
                << " defined in the dynamicMeshDict dictionary"
                << abort(FatalError);
        }

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        forAll(eFaces, edgeI)
        {
            deltaHf[eFaces[edgeI]] *= 2.0;
        }
    }

    controlPoints() += facesDisplacementDir()*deltaHf;

    displacement = pointDisplacement();

    velocityMotionSolver& vMotion =
        refCast<velocityMotionSolver>
        (
            const_cast<motionSolver&>(motion())
        );

    pointVectorField& pointMotionU = vMotion.pointMotionU();
    pointMotionU.primitiveFieldRef() = Zero;

    fixedValuePointPatchVectorField& fsPatchPointMeshU =
        refCast<fixedValuePointPatchVectorField>
        (
            const_cast<pointPatchVectorField&>
            (
                pointMotionU.boundaryField()[fsPatchIndex()]
            )
        );

    fsPatchPointMeshU == displacement/mesh().time().deltaTValue();

    dynamicMotionSolverFvMesh::update();
}


void Foam::interfaceTrackingFvMesh::updateSurfaceFlux()
{
    Phis() = fac::interpolate(Us()) & aMesh().Le();
}


void Foam::interfaceTrackingFvMesh::updateSurfactantConcentration()
{
    if (!pureFreeSurface())
    {
        Info<< "Correct surfactant concentration" << endl << flush;

        updateSurfaceFlux();

        // Crate and solve the surfactanta transport equation
        faScalarMatrix CsEqn
        (
            fam::ddt(surfactantConcentration())
          + fam::div(Phis(), surfactantConcentration())
          - fam::laplacian
            (
                surfactant().diffusion(),
                surfactantConcentration()
            )
        );

        if (surfactant().soluble())
        {
            #include "solveBulkSurfactant.H"

            areaScalarField Cb
            (
                IOobject
                (
                    "Cb",
                    aMesh().time().timeName(),
                    aMesh().thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh(),
                dimensionedScalar(dimMoles/dimVolume, Zero),
                faPatchFieldBase::zeroGradientType()
            );

            Cb.ref().field() =
                bulkSurfactantConcentration().boundaryField()[fsPatchIndex()];
            Cb.correctBoundaryConditions();

            CsEqn +=
                fam::Sp
                (
                    surfactant().adsorptionCoeff()*Cb
                  + surfactant().adsorptionCoeff()
                   *surfactant().desorptionCoeff(),
                    surfactantConcentration()
                )
              - surfactant().adsorptionCoeff()
               *Cb*surfactant().saturatedConc();
        }

        CsEqn.solve();

        // Info<< "Correct surface tension" << endl;

        surfaceTension() =
            sigma() + surfactant().dSigma(surfactantConcentration())/rho_;

        if (neg(min(surfaceTension().internalField().field())))
        {
            FatalErrorInFunction
                << "Surface tension is negative"
                << abort(FatalError);
        }
    }
}


Foam::vector Foam::interfaceTrackingFvMesh::totalPressureForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& P = p().boundaryField()[fsPatchIndex()];

    vectorField pressureForces(S*P*n);

    return gSum(pressureForces);
}


Foam::vector Foam::interfaceTrackingFvMesh::totalViscousForce() const
{
    const auto& turbulence =
        mesh().lookupObject<turbulenceModel>("turbulenceProperties");

    scalarField nu(turbulence.nuEff(fsPatchIndex()));

    // const singlePhaseTransportModel& properties =
    //     mesh().lookupObject<singlePhaseTransportModel>
    //     (
    //         "transportProperties"
    //     );

    // dimensionedScalar nu("nu", properties);

    const scalarField& S = aMesh().S();
    const vectorField& n = aMesh().faceAreaNormals().internalField();

    vectorField nGradU
    (
        U().boundaryField()[fsPatchIndex()].snGrad()
    );

    vectorField viscousForces
    (
      - nu*S
       *(
            nGradU
          + (fac::grad(Us())().internalField()&n)
          - (n*fac::div(Us())().internalField())
        )
    );

    return gSum(viscousForces);
}


Foam::vector Foam::interfaceTrackingFvMesh::totalSurfaceTensionForce() const
{
    const scalarField& S = aMesh().S();

    const vectorField& n = aMesh().faceAreaNormals().internalField();

    const scalarField& K = aMesh().faceCurvatures().internalField();

    vectorField surfTensionForces(n.size(), Zero);

    if (pureFreeSurface())
    {
        surfTensionForces =
            S*sigma().value()
           *fac::edgeIntegrate
            (
                aMesh().Le()*aMesh().edgeLengthCorrection()
            )().internalField();
    }
    else
    {
        surfTensionForces = surfaceTension().internalField().field()*K*S*n;
    }

    return gSum(surfTensionForces);
}


Foam::scalar Foam::interfaceTrackingFvMesh::maxCourantNumber()
{
    scalar CoNum = 0;

    if (pureFreeSurface())
    {
        const scalarField& dE = aMesh().lPN();

        CoNum = gMax
        (
            mesh().time().deltaTValue()/
            sqrt
            (
                Foam::pow(dE, 3.0)/2.0/M_PI/(sigma().value() + SMALL)
            )
        );
    }
    else
    {
        scalarField sigmaE
        (
            linearEdgeInterpolate(surfaceTension())().internalField().field()
          + SMALL
        );

        const scalarField& dE = aMesh().lPN();

        CoNum = gMax
        (
            mesh().time().deltaTValue()/
            sqrt
            (
                Foam::pow(dE, 3.0)/2.0/M_PI/sigmaE
            )
        );
    }

    return CoNum;
}


void Foam::interfaceTrackingFvMesh::updateProperties()
{
    const singlePhaseTransportModel& properties =
        mesh().lookupObject<singlePhaseTransportModel>
        (
            "transportProperties"
        );

    rho_ = dimensionedScalar("rho", properties);

    sigma0_ = dimensionedScalar("sigma", properties)/rho_;
}


void Foam::interfaceTrackingFvMesh::correctPointDisplacement
(
    const scalarField& sweptVolCorr,
    vectorField& displacement
)
{
    const labelListList& pFaces =
        aMesh().patch().pointFaces();

    const faceList& faces =
        aMesh().patch().localFaces();

    const pointField& points =
        aMesh().patch().localPoints();

    for (const word& patchName : fixedFreeSurfacePatches_)
    {
        label fixedPatchID = aMesh().boundary().findPatchID(patchName);

        const labelList& pLabels =
            aMesh().boundary()[fixedPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[fixedPatchID].edgeFaces();

        labelHashSet pointSet;

        forAll(eFaces, edgeI)
        {
            label curFace = eFaces[edgeI];

            const labelList& curPoints = faces[curFace];

            forAll(curPoints, pointI)
            {
                label curPoint = curPoints[pointI];
                label index = pLabels.find(curPoint);

                if (index == -1)
                {
                    pointSet.insert(curPoint);
                }
            }
        }

        labelList corrPoints = pointSet.toc();

        labelListList corrPointFaces(corrPoints.size());

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            labelHashSet faceSet;

            forAll(pFaces[curPoint], faceI)
            {
                label curFace = pFaces[curPoint][faceI];

                label index = eFaces.find(curFace);

                if (index != -1)
                {
                    faceSet.insert(curFace);
                }
            }

            corrPointFaces[pointI] = faceSet.toc();
        }

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            scalar curDisp = 0;

            const labelList& curPointFaces = corrPointFaces[pointI];

            forAll(curPointFaces, faceI)
            {
                const face& curFace = faces[curPointFaces[faceI]];

                label ptInFace = curFace.which(curPoint);
                label next = curFace.nextLabel(ptInFace);
                label prev = curFace.prevLabel(ptInFace);

                vector a = points[next] - points[curPoint];
                vector b = points[prev] - points[curPoint];
                const vector& c = pointsDisplacementDir()[curPoint];

                curDisp += 2*sweptVolCorr[curPointFaces[faceI]]/((a^b)&c);
            }

            curDisp /= curPointFaces.size();

            displacement[curPoint] =
                curDisp*pointsDisplacementDir()[curPoint];
        }
    }


    for (const word& patchName : nonReflectingFreeSurfacePatches_)
    {
        label nonReflectingPatchID =
            aMesh().boundary().findPatchID(patchName);

        const labelList& pLabels =
            aMesh().boundary()[nonReflectingPatchID].pointLabels();

        const labelList& eFaces =
            aMesh().boundary()[nonReflectingPatchID].edgeFaces();

        labelList corrPoints = pLabels;

        labelListList corrPointFaces(corrPoints.size());

        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            labelHashSet faceSet;

            forAll(pFaces[curPoint], faceI)
            {
                label curFace = pFaces[curPoint][faceI];

                label index = eFaces.find(curFace);

                if (index != -1)
                {
                    faceSet.insert(curFace);
                }
            }

            corrPointFaces[pointI] = faceSet.toc();
        }


        forAll(corrPoints, pointI)
        {
            label curPoint = corrPoints[pointI];

            scalar curDisp = 0;

            const labelList& curPointFaces = corrPointFaces[pointI];

            forAll(curPointFaces, faceI)
            {
                const face& curFace = faces[curPointFaces[faceI]];

                label ptInFace = curFace.which(curPoint);
                label next = curFace.nextLabel(ptInFace);
                label prev = curFace.prevLabel(ptInFace);

                label p0 = -1;
                label p1 = -1;
                label p2 = -1;

                if (corrPoints.find(next) == -1)
                {
                    p0 = curPoint;
                    p1 = next;
                    p2 = curFace.nextLabel(curFace.which(next));
                }
                else
                {
                    p0 = curFace.prevLabel(curFace.which(prev));
                    p1 = prev;
                    p2 = curPoint;
                }

                vector a0 = points[p1] - points[p0];
                vector b0 = points[p2] - points[p1];
                vector c0 = displacement[p1];

                scalar V0 = mag((a0^b0)&c0)/2;

                scalar DV = sweptVolCorr[curPointFaces[faceI]] - V0;

                if (corrPoints.find(prev) != -1)
                {
                    vector a = points[curPoint] - points[prev];
                    vector b =
                        (points[next] + displacement[next])
                      - points[curPoint];
                    const vector& c = pointsDisplacementDir()[curPoint];

                    curDisp += 2*DV/((a^b)&c);
                }
                else
                {
                    vector a = points[curPoint]
                      - (points[prev] + displacement[prev]);
                    vector b = points[next] - points[curPoint];
                    const vector& c = pointsDisplacementDir()[curPoint];

                    curDisp += 2*DV/((a^b)&c);
                }
            }

            curDisp /= curPointFaces.size();

            displacement[curPoint] =
                curDisp*pointsDisplacementDir()[curPoint];
        }
    }
}


void Foam::interfaceTrackingFvMesh::correctContactLinePointNormals()
{
    // Correct normals for contact line points
    // according to specified contact angle

    vectorField& N =
        const_cast<vectorField&>
        (
            aMesh().pointAreaNormals()
        );

    if (contactAnglePtr_ && correctContactLineNormals())
    {
        Info<< "Correcting contact line normals" << endl;

        vectorField oldPoints(aMesh().nPoints(), Zero);

        const labelList& meshPoints = aMesh().patch().meshPoints();

        forAll(oldPoints, ptI)
        {
            oldPoints[ptI] =
                mesh().oldPoints()[meshPoints[ptI]];
        }

// #       include "createTangentField.H"
        areaVectorField tangent
        (
            IOobject
            (
                "tangent",
                aMesh().time().timeName(),
                aMesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh(),
            dimensioned<vector>(dimless, Zero)
        );

        if (Pstream::parRun())
        {
            const labelListList& edgeFaces = aMesh().patch().edgeFaces();
            const labelListList& pointEdges = aMesh().patch().pointEdges();
            const labelListList& pointFaces = aMesh().patch().pointFaces();
            const edgeList& edges = aMesh().edges();

            forAll(aMesh().boundary(), patchI)
            {
                if
                (
                    aMesh().boundary()[patchI].type()
                 == processorFaPatch::typeName
                )
                {
                    const processorFaPatch& procPatch =
                        refCast<const processorFaPatch>
                        (
                            aMesh().boundary()[patchI]
                        );

                    const labelList& patchPointLabels =
                        procPatch.pointLabels();

                    forAll(patchPointLabels, pointI)
                    {
                        label curPoint = patchPointLabels[pointI];

                        // Check if processor point is boundary point

                        label patchID = -1;
                        label edgeID = -1;

                        const labelList& curPointEdges = pointEdges[curPoint];

                        forAll(curPointEdges, edgeI)
                        {
                            label curEdge = curPointEdges[edgeI];

                            if (edgeFaces[curEdge].size() == 1)
                            {
                                forAll(aMesh().boundary(), pI)
                                {
                                    const labelList& curEdges =
                                        aMesh().boundary()[pI];

                                    label index = curEdges.find(curEdge);

                                    if (index != -1)
                                    {
                                        if
                                        (
                                            aMesh().boundary()[pI].type()
                                         != processorFaPatch::typeName
                                        )
                                        {
                                            patchID = pI;
                                            edgeID = index;
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        if (patchID != -1)
                        {
                            label curEdge =
                                aMesh().boundary()[patchID].start() + edgeID;

                            vector t = edges[curEdge].vec(oldPoints);
                            t /= mag(t) + SMALL;

                            const labelList& curPointFaces =
                                pointFaces[curPoint];

                            forAll(curPointFaces, fI)
                            {
                                tangent.ref().field()[curPointFaces[fI]] = t;
                            }
                        }
                    }
                }
            }

            tangent.correctBoundaryConditions();
        }

        forAll(aMesh().boundary(), patchI)
        {
            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    mesh().boundary()[ngbPolyPatchID].type()
                 == wallFvPatch::typeName
                )
                {
                    const scalar rotAngle = degToRad
                    (
                        gAverage
                        (
                            90
                          - contactAnglePtr_->boundaryField()[patchI]
                        )
                    );

                    vectorField ngbN
                    (
                        aMesh().boundary()[patchI].ngbPolyPatchPointNormals()
                    );

                    const labelList& patchPoints =
                        aMesh().boundary()[patchI].pointLabels();

                    vectorField pN(N, patchPoints);

                    vectorField rotationAxis(ngbN^pN);
                    rotationAxis /= mag(rotationAxis) + SMALL;


                    // Calc rotation axis using edge vectors

                    const edgeList& edges = aMesh().edges();

                    const labelListList& pointEdges =
                        aMesh().boundary()[patchI].pointEdges();

                    forAll(pointEdges, pointI)
                    {
                        vector rotAx = Zero;

                        forAll(pointEdges[pointI], eI)
                        {
                            label curEdge =
                                aMesh().boundary()[patchI].start()
                              + pointEdges[pointI][eI];

                            vector e = edges[curEdge].vec(oldPoints);

                            e *= (e&rotationAxis[pointI])
                               /mag(e&rotationAxis[pointI]);

                            e /= mag(e) + SMALL;

                            rotAx += e;
                        }

                        if (pointEdges[pointI].size() == 1)
                        {
                            label curPoint = patchPoints[pointI];

                            const labelListList& ptEdges =
                                aMesh().patch().pointEdges();
                            const labelList& curPointEdges =
                                ptEdges[curPoint];

                            label procPatchID = -1;
                            label edgeID = -1;

                            const labelListList& edgeFaces =
                                aMesh().patch().edgeFaces();

                            forAll(curPointEdges, edgeI)
                            {
                                label curEdge = curPointEdges[edgeI];

                                if (edgeFaces[curEdge].size() == 1)
                                {
                                    forAll(aMesh().boundary(), pI)
                                    {
                                        const labelList& curEdges =
                                            aMesh().boundary()[pI];

                                        label index =
                                            curEdges.find(curEdge);

                                        if (index != -1)
                                        {
                                            if
                                            (
                                                aMesh().boundary()[pI].type()
                                                == processorFaPatch::typeName
                                            )
                                            {
                                                procPatchID = pI;
                                                edgeID = index;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }

                            if (procPatchID != -1)
                            {
                                vector t =
                                    tangent.boundaryField()[procPatchID]
                                   .patchNeighbourField()()[edgeID];

                                t *= (t&rotationAxis[pointI])
                                    /(mag(t&rotationAxis[pointI]) + SMALL);

                                t /= mag(t) + SMALL;

                                rotAx += t;
                            }
                        }

                        rotationAxis[pointI] = rotAx/(mag(rotAx) + SMALL);
                    }

                    // Rodrigues' rotation formula
                    ngbN = ngbN*cos(rotAngle)
                      + rotationAxis*(rotationAxis & ngbN)*(1 - cos(rotAngle))
                      + (rotationAxis^ngbN)*sin(rotAngle);

                    // Info<< aMesh().boundary()[patchI].name() << endl;
                    forAll(patchPoints, pointI)
                    {
                        N[patchPoints[pointI]] -=
                            ngbN[pointI]*(ngbN[pointI]&N[patchPoints[pointI]]);

                        N[patchPoints[pointI]] /=
                            mag(N[patchPoints[pointI]]) + SMALL;

                        // Info<< N[patchPoints[pointI]] << endl;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceTrackingFvMesh::interfaceTrackingFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicMotionSolverFvMesh(io, doInit),
    aMeshPtr_(nullptr),
    fsPatchIndex_(-1),
    fixedFreeSurfacePatches_(),
    nonReflectingFreeSurfacePatches_(),
    pointNormalsCorrectionPatches_(),
    normalMotionDir_(false),
    motionDir_(Zero),
    smoothing_(false),
    pureFreeSurface_(true),
    rigidFreeSurface_(false),
    correctContactLineNormals_(false),
    sigma0_("zero", dimForce/dimLength/dimDensity, Zero),
    rho_("one", dimDensity, 1.0),
    timeIndex_(-1),
    UsPtr_(nullptr),
    controlPointsPtr_(nullptr),
    motionPointsMaskPtr_(nullptr),
    pointsDisplacementDirPtr_(nullptr),
    facesDisplacementDirPtr_(nullptr),
    fsNetPhiPtr_(nullptr),
    phisPtr_(nullptr),
    surfactConcPtr_(nullptr),
    bulkSurfactConcPtr_(nullptr),
    surfaceTensionPtr_(nullptr),
    surfactantPtr_(nullptr),
    contactAnglePtr_(nullptr)
{
    if (doInit)
    {
        init(false);    // do not initialise lower levels
    }
}

/*
Foam::interfaceTrackingFvMesh::interfaceTrackingFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    dynamicMotionSolverFvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    aMeshPtr_(new faMesh(*this)),
    fsPatchIndex_(-1),
    fixedFreeSurfacePatches_(),
    nonReflectingFreeSurfacePatches_(),
    pointNormalsCorrectionPatches_(),
    normalMotionDir_(false),
    motionDir_(Zero),
    smoothing_(false),
    pureFreeSurface_(true),
    sigma0_("zero", dimForce/dimLength/dimDensity, Zero),
    rho_("one", dimDensity, 1.0),
    timeIndex_(-1),
    UsPtr_(nullptr),
    controlPointsPtr_(nullptr),
    motionPointsMaskPtr_(nullptr),
    pointsDisplacementDirPtr_(nullptr),
    facesDisplacementDirPtr_(nullptr),
    fsNetPhiPtr_(nullptr),
    phisPtr_(nullptr),
    surfactConcPtr_(nullptr),
    bulkSurfactConcPtr_(nullptr),
    surfaceTensionPtr_(nullptr),
    surfactantPtr_(nullptr),
    contactAnglePtr_(nullptr)
{}
*/

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceTrackingFvMesh::~interfaceTrackingFvMesh()
{
    deleteDemandDrivenData(UsPtr_);
    deleteDemandDrivenData(controlPointsPtr_);
    deleteDemandDrivenData(motionPointsMaskPtr_);
    deleteDemandDrivenData(pointsDisplacementDirPtr_);
    deleteDemandDrivenData(facesDisplacementDirPtr_);
    deleteDemandDrivenData(fsNetPhiPtr_);
    deleteDemandDrivenData(phisPtr_);
    deleteDemandDrivenData(surfactConcPtr_);
    deleteDemandDrivenData(bulkSurfactConcPtr_);
    deleteDemandDrivenData(surfaceTensionPtr_);
    deleteDemandDrivenData(surfactantPtr_);
    deleteDemandDrivenData(contactAnglePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::interfaceTrackingFvMesh::init(const bool doInit)
{
    if (doInit)
    {
        dynamicMotionSolverFvMesh::init(doInit);
    }

    aMeshPtr_.reset(new faMesh(*this));

    // Set motion-based data
    fixedFreeSurfacePatches_ =
        motion().get<wordList>("fixedFreeSurfacePatches");

    pointNormalsCorrectionPatches_ =
        motion().get<wordList>("pointNormalsCorrectionPatches");

    normalMotionDir_ = motion().get<bool>("normalMotionDir");
    smoothing_ = motion().getOrDefault("smoothing", false);
    pureFreeSurface_ = motion().getOrDefault("pureFreeSurface", true);

    initializeData();

    return true;
}


Foam::areaVectorField& Foam::interfaceTrackingFvMesh::Us()
{
    if (!UsPtr_)
    {
        makeUs();
    }

    return *UsPtr_;
}


const Foam::areaVectorField& Foam::interfaceTrackingFvMesh::Us() const
{
    if (!UsPtr_)
    {
        makeUs();
    }

    return *UsPtr_;
}


Foam::areaScalarField& Foam::interfaceTrackingFvMesh::fsNetPhi()
{
    if (!fsNetPhiPtr_)
    {
        makeFsNetPhi();
    }

    return *fsNetPhiPtr_;
}


const Foam::areaScalarField& Foam::interfaceTrackingFvMesh::fsNetPhi() const
{
    if (!fsNetPhiPtr_)
    {
        makeFsNetPhi();
    }

    return *fsNetPhiPtr_;
}


void Foam::interfaceTrackingFvMesh::correctUsBoundaryConditions()
{
    forAll(Us().boundaryField(), patchI)
    {
        if
        (
            Us().boundaryField()[patchI].type()
         == calculatedFaPatchVectorField::typeName
        )
        {
            vectorField& pUs = Us().boundaryFieldRef()[patchI];

            pUs = Us().boundaryField()[patchI].patchInternalField();

            label ngbPolyPatchID =
                aMesh().boundary()[patchI].ngbPolyPatchIndex();

            if (ngbPolyPatchID != -1)
            {
                if
                (
                    (
                        U().boundaryField()[ngbPolyPatchID].type()
                     == slipFvPatchVectorField::typeName
                    )
                 ||
                    (
                        U().boundaryField()[ngbPolyPatchID].type()
                     == symmetryFvPatchVectorField::typeName
                    )
                )
                {
                    vectorField N
                    (
                        aMesh().boundary()[patchI].ngbPolyPatchFaceNormals()
                    );

                    pUs -= N*(N&pUs);
                }
            }
        }
    }

    Us().correctBoundaryConditions();
}


void Foam::interfaceTrackingFvMesh::updateUs()
{
    // Info<< "Update Us" << endl;

    Us().ref().field() = U().boundaryField()[fsPatchIndex()];

    // // Correct normal component of free-surface velocity
    // const vectorField& nA = aMesh().faceAreaNormals().internalField();
    // vectorField UnFs = nA*phi().boundaryField()[fsPatchIndex()]
    //    /mesh().boundary()[fsPatchIndex()].magSf();
    // Us().ref().field() += UnFs - nA*(nA&Us().internalField());

    correctUsBoundaryConditions();
}


const Foam::volVectorField& Foam::interfaceTrackingFvMesh::U() const
{
    return *getObjectPtr<const volVectorField>("U");
}


const Foam::volScalarField& Foam::interfaceTrackingFvMesh::p() const
{
    return *getObjectPtr<const volScalarField>("p");
}


const Foam::surfaceScalarField& Foam::interfaceTrackingFvMesh::phi() const
{
    return *getObjectPtr<const surfaceScalarField>("phi");
}


Foam::tmp<Foam::vectorField>
Foam::interfaceTrackingFvMesh::freeSurfaceSnGradU()
{
    auto tSnGradU = tmp<vectorField>::New(aMesh().nFaces(), Zero);
    auto& SnGradU = tSnGradU.ref();

    const vectorField& nA = aMesh().faceAreaNormals().internalField();

    areaScalarField divUs
    (
        fac::div(Us())
      - aMesh().faceCurvatures()*(aMesh().faceAreaNormals()&Us())
    );

    areaTensorField gradUs(fac::grad(Us()));

    // Remove component of gradient normal to surface (area)
    const areaVectorField& n = aMesh().faceAreaNormals();
    gradUs -= n*(n & gradUs);
    gradUs.correctBoundaryConditions();

    const turbulenceModel& turbulence =
        mesh().lookupObject<turbulenceModel>("turbulenceProperties");

    scalarField nu(turbulence.nuEff(fsPatchIndex()));

    vectorField tangentialSurfaceTensionForce(nA.size(), Zero);

    if (!pureFreeSurface() && max(nu) > SMALL)
    {
        tangentialSurfaceTensionForce =
            surfaceTensionGrad()().internalField();
    }

    SnGradU =
        tangentialSurfaceTensionForce/(nu + SMALL)
      - nA*divUs.internalField()
      - (gradUs.internalField()&nA);

    return tSnGradU;
}


Foam::tmp<Foam::scalarField>
Foam::interfaceTrackingFvMesh::freeSurfaceSnGradUn()
{
    auto tSnGradUn = tmp<scalarField>::New(aMesh().nFaces(), Zero);
    auto& SnGradUn = tSnGradUn.ref();

    areaScalarField divUs
    (
        fac::div(Us())
      - aMesh().faceCurvatures()*(aMesh().faceAreaNormals()&Us())
    );

    SnGradUn = -divUs.internalField();

    return tSnGradUn;
}


Foam::tmp<scalarField>
Foam::interfaceTrackingFvMesh::freeSurfacePressureJump()
{
    auto tPressureJump = tmp<scalarField>::New(aMesh().nFaces(), Zero);
    auto& pressureJump = tPressureJump.ref();

    const scalarField& K = aMesh().faceCurvatures().internalField();

    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(mesh().time());

    const turbulenceModel& turbulence =
        mesh().lookupObject<turbulenceModel>("turbulenceProperties");

    scalarField nu(turbulence.nuEff(fsPatchIndex()));

    pressureJump =
      - (g.value() & mesh().Cf().boundaryField()[fsPatchIndex()])
      + 2.0*nu*freeSurfaceSnGradUn();

    if (pureFreeSurface())
    {
        pressureJump -= sigma().value()*K;
    }
    else
    {
        pressureJump -= surfaceTension().internalField()*K;
    }

    return tPressureJump;
}


Foam::vectorField& Foam::interfaceTrackingFvMesh::controlPoints()
{
    if (!controlPointsPtr_)
    {
        makeControlPoints();
    }

    return *controlPointsPtr_;
}


Foam::labelList& Foam::interfaceTrackingFvMesh::motionPointsMask()
{
    if (!motionPointsMaskPtr_)
    {
        makeMotionPointsMask();
    }

    return *motionPointsMaskPtr_;
}


Foam::vectorField& Foam::interfaceTrackingFvMesh::pointsDisplacementDir()
{
    if (!pointsDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *pointsDisplacementDirPtr_;
}


Foam::vectorField& Foam::interfaceTrackingFvMesh::facesDisplacementDir()
{
    if (!facesDisplacementDirPtr_)
    {
        makeDirections();
    }

    return *facesDisplacementDirPtr_;
}


Foam::edgeScalarField& Foam::interfaceTrackingFvMesh::Phis()
{
    if (!phisPtr_)
    {
        makePhis();
    }

    return *phisPtr_;
}

Foam::areaScalarField&
Foam::interfaceTrackingFvMesh::surfactantConcentration()
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }

    return *surfactConcPtr_;
}


const Foam::areaScalarField&
Foam::interfaceTrackingFvMesh::surfactantConcentration() const
{
    if (!surfactConcPtr_)
    {
        makeSurfactConc();
    }

    return *surfactConcPtr_;
}


Foam::volScalarField&
Foam::interfaceTrackingFvMesh::bulkSurfactantConcentration()
{
    if (!bulkSurfactConcPtr_)
    {
        makeBulkSurfactConc();
    }

    return *bulkSurfactConcPtr_;
}


const Foam::volScalarField&
Foam::interfaceTrackingFvMesh::bulkSurfactantConcentration() const
{
    if (!bulkSurfactConcPtr_)
    {
        makeBulkSurfactConc();
    }

    return *bulkSurfactConcPtr_;
}


Foam::areaScalarField&
Foam::interfaceTrackingFvMesh::surfaceTension()
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }

    return *surfaceTensionPtr_;
}


const Foam::areaScalarField&
Foam::interfaceTrackingFvMesh::surfaceTension() const
{
    if (!surfaceTensionPtr_)
    {
        makeSurfaceTension();
    }

    return *surfaceTensionPtr_;
}


const Foam::surfactantProperties&
Foam::interfaceTrackingFvMesh::surfactant() const
{
    if (!surfactantPtr_)
    {
        makeSurfactant();
    }

    return *surfactantPtr_;
}


Foam::tmp<Foam::areaVectorField>
Foam::interfaceTrackingFvMesh::surfaceTensionGrad()
{
    auto tgrad = tmp<areaVectorField>::New
    (
        IOobject
        (
            "surfaceTensionGrad",
            aMesh().time().timeName(),
            aMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fac::grad(surfaceTension())
        // (-fac::grad(surfactantConcentration()/rho_)*
        // surfactant().surfactR()*surfactant().surfactT()/
        // (1.0 - surfactantConcentration()/
        // surfactant().surfactSaturatedConc()))()
    );
    auto& grad = tgrad.ref();

    // Remove component of gradient normal to surface (area)
    const areaVectorField& n = aMesh().faceAreaNormals();
    grad -= n*(n & grad);
    grad.correctBoundaryConditions();

    return tgrad;
}


bool Foam::interfaceTrackingFvMesh::update()
{
    if (timeIndex_ != mesh().time().timeIndex())
    {
        if (smoothing_ && !rigidFreeSurface_)
        {
            smoothFreeSurfaceMesh();
            clearControlPoints();
        }

        updateDisplacementDirections();

        updateProperties();

        Info<< "Maximal capillary Courant number: "
            << maxCourantNumber() << endl;

        const scalarField& K = aMesh().faceCurvatures().internalField();

        Info<< "Free surface curvature: min = " << gMin(K)
            << ", max = " << gMax(K) << ", average = " << gAverage(K) << nl;

        timeIndex_ = mesh().time().timeIndex();
    }

    if (!rigidFreeSurface_)
    {
        // This is currently relaltive flux
        scalarField sweptVolCorr =
            phi().boundaryField()[fsPatchIndex()];

        // Info<< "Free surface flux: sum local = "
        //     << gSum(mag(sweptVolCorr))
        //     << ", global = " << gSum(sweptVolCorr) << endl;

        // if (mesh().moving())
        // {
        //     sweptVolCorr -=
        //         fvc::meshPhi(U())().boundaryField()[fsPatchIndex()];
        // }

        Info<< "Free surface continuity error : sum local = "
            << gSum(mag(sweptVolCorr)) << ", global = " << gSum(sweptVolCorr)
            << endl;

        // For postprocessing
        fsNetPhi().ref().field() = sweptVolCorr;

        word ddtScheme
        (
            mesh().ddtScheme
            (
                "ddt(" + U().name() + ')'
            )
        );

        if
        (
            ddtScheme
         == fv::CrankNicolsonDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= (1.0/2.0)*mesh().time().deltaTValue();
        }
        else if
        (
            ddtScheme
         == fv::EulerDdtScheme<vector>::typeName
        )
        {
            sweptVolCorr *= mesh().time().deltaTValue();
        }
        else if
        (
            ddtScheme
         == fv::backwardDdtScheme<vector>::typeName
        )
        {
            if (mesh().time().timeIndex() == 1)
            {
                sweptVolCorr *= mesh().time().deltaTValue();
            }
            else
            {
                sweptVolCorr *= (2.0/3.0)*mesh().time().deltaTValue();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported temporal differencing scheme : "
                << ddtScheme << nl
                << abort(FatalError);
        }

        const scalarField& Sf = aMesh().S();
        const vectorField& Nf = aMesh().faceAreaNormals().internalField();

        scalarField deltaHf
        (
            sweptVolCorr/(Sf*(Nf & facesDisplacementDir()))
        );

        for (const word& patchName : fixedFreeSurfacePatches_)
        {
            label fixedPatchID =
                aMesh().boundary().findPatchID(patchName);

            if (fixedPatchID == -1)
            {
                FatalErrorInFunction
                    << "Wrong faPatch name '" << patchName
                    << "'in the fixedFreeSurfacePatches list"
                    << " defined in dynamicMeshDict dictionary"
                    << abort(FatalError);
            }

            const labelList& eFaces =
                aMesh().boundary()[fixedPatchID].edgeFaces();

            forAll(eFaces, edgeI)
            {
                deltaHf[eFaces[edgeI]] *= 2.0;
            }
        }

        controlPoints() += facesDisplacementDir()*deltaHf;

        pointField displacement(pointDisplacement());
        correctPointDisplacement(sweptVolCorr, displacement);

        velocityMotionSolver& vMotion =
            refCast<velocityMotionSolver>
            (
                const_cast<motionSolver&>(motion())
            );

        pointVectorField& pointMotionU = vMotion.pointMotionU();
        pointMotionU.primitiveFieldRef() = Zero;

        fixedValuePointPatchVectorField& fsPatchPointMeshU =
            refCast<fixedValuePointPatchVectorField>
            (
                const_cast<pointPatchVectorField&>
                (
                    pointMotionU.boundaryField()[fsPatchIndex()]
                )
            );

        fsPatchPointMeshU == displacement/mesh().time().deltaTValue();

        dynamicMotionSolverFvMesh::update();

        correctContactLinePointNormals();
    }
    else
    {
        vectorField displacement
        (
            mesh().boundaryMesh()[fsPatchIndex()].nPoints(),
            Zero
        );

        velocityMotionSolver& vMotion =
            refCast<velocityMotionSolver>
            (
                const_cast<motionSolver&>(motion())
            );

        pointVectorField& pointMotionU = vMotion.pointMotionU();
        pointMotionU.primitiveFieldRef() = Zero;

        fixedValuePointPatchVectorField& fsPatchPointMeshU =
            refCast<fixedValuePointPatchVectorField>
            (
                const_cast<pointPatchVectorField&>
                (
                    pointMotionU.boundaryField()[fsPatchIndex()]
                )
            );

        fsPatchPointMeshU == displacement/mesh().time().deltaTValue();

        dynamicMotionSolverFvMesh::update();
    }

    updateUs();

    updateSurfactantConcentration();

    return true;
}


void Foam::interfaceTrackingFvMesh::writeVTK() const
{
    vtk::uindirectPatchWriter writer
    (
        aMesh().patch(),
        vtk::formatType::LEGACY_ASCII,
        mesh().time().timePath()/"freeSurface",
        false // serial only
    );
    writer.writeGeometry();
}


void Foam::interfaceTrackingFvMesh::writeVTKControlPoints()
{
    // Write control points into VTK
    OFstream os
    (
        mesh().time().timePath()/"freeSurfaceControlPoints.vtk"
    );

    Info<< "Writing free surface control points to " << os.name() << nl;

    os  << "# vtk DataFile Version 2.0" << nl
        << "freeSurfaceControlPoints" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    const label nPoints = controlPoints().size();

    os  << "POINTS " << nPoints << " float" << nl;
    for (const point& p : controlPoints())
    {
        os  << float(p.x()) << ' '
            << float(p.y()) << ' '
            << float(p.z()) << nl;
    }

    os  << "VERTICES " << nPoints << ' ' << 2*nPoints << nl;
    for (label id = 0; id < nPoints; ++id)
    {
        os  << 1 << ' ' << id << nl;
    }
}


// ************************************************************************* //

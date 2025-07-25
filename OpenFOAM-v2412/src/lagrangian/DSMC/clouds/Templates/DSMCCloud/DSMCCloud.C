/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "DSMCCloud.H"
#include "BinaryCollisionModel.H"
#include "WallInteractionModel.H"
#include "InflowBoundaryModel.H"
#include "constants.H"
#include "zeroGradientFvPatchFields.H"
#include "polyMeshTetDecomposition.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::buildConstProps()
{
    Info<< nl << "Constructing constant properties for" << endl;
    constProps_.setSize(typeIdList_.size());

    dictionary moleculeProperties
    (
        particleProperties_.subDict("moleculeProperties")
    );

    forAll(typeIdList_, i)
    {
        const word& id = typeIdList_[i];

        Info<< "    " << id << endl;

        const dictionary& molDict(moleculeProperties.subDict(id));

        constProps_[i] = typename ParcelType::constantProperties(molDict);
    }
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::buildCellOccupancy()
{
    for (auto& list : cellOccupancy_)
    {
        list.clear();
    }

    for (ParcelType& p : *this)
    {
        cellOccupancy_[p.cell()].append(&p);
    }
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::initialise
(
    const IOdictionary& dsmcInitialiseDict
)
{
    Info<< nl << "Initialising particles" << endl;

    const scalar temperature
    (
        dsmcInitialiseDict.get<scalar>("temperature")
    );

    const vector velocity(dsmcInitialiseDict.get<vector>("velocity"));

    const dictionary& numberDensitiesDict
    (
        dsmcInitialiseDict.subDict("numberDensities")
    );

    List<word> molecules(numberDensitiesDict.toc());

    Field<scalar> numberDensities(molecules.size());

    forAll(molecules, i)
    {
        numberDensities[i] = numberDensitiesDict.get<scalar>(molecules[i]);
    }

    numberDensities /= nParticle_;

    forAll(mesh_.cells(), celli)
    {
        List<tetIndices> cellTets = polyMeshTetDecomposition::cellTetIndices
        (
            mesh_,
            celli
        );

        forAll(cellTets, tetI)
        {
            const tetIndices& cellTetIs = cellTets[tetI];
            tetPointRef tet = cellTetIs.tet(mesh_);
            scalar tetVolume = tet.mag();

            forAll(molecules, i)
            {
                const word& moleculeName(molecules[i]);

                label typeId = typeIdList_.find(moleculeName);

                if (typeId == -1)
                {
                    FatalErrorInFunction
                        << "typeId " << moleculeName << "not defined." << nl
                        << abort(FatalError);
                }

                const typename ParcelType::constantProperties& cP =
                    constProps(typeId);

                scalar numberDensity = numberDensities[i];

                // Calculate the number of particles required
                scalar particlesRequired = numberDensity*tetVolume;

                // Only integer numbers of particles can be inserted
                label nParticlesToInsert = label(particlesRequired);

                // Add another particle with a probability proportional to the
                // remainder of taking the integer part of particlesRequired
                if
                (
                    (particlesRequired - nParticlesToInsert)
                  > rndGen_.sample01<scalar>()
                )
                {
                    nParticlesToInsert++;
                }

                for (label pI = 0; pI < nParticlesToInsert; pI++)
                {
                    point p = tet.randomPoint(rndGen_);

                    vector U = equipartitionLinearVelocity
                    (
                        temperature,
                        cP.mass()
                    );

                    scalar Ei = equipartitionInternalEnergy
                    (
                        temperature,
                        cP.internalDegreesOfFreedom()
                    );

                    U += velocity;

                    addNewParcel(p, celli, U, Ei, typeId);
                }
            }
        }
    }

    // Initialise the sigmaTcRMax_ field to the product of the cross section of
    // the most abundant species and the most probable thermal speed (Bird,
    // p222-223)

    label mostAbundantType(findMax(numberDensities));

    const typename ParcelType::constantProperties& cP = constProps
    (
        mostAbundantType
    );

    sigmaTcRMax_.primitiveFieldRef() = cP.sigmaT()*maxwellianMostProbableSpeed
    (
        temperature,
        cP.mass()
    );

    sigmaTcRMax_.correctBoundaryConditions();
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::collisions()
{
    if (!binaryCollision().active())
    {
        return;
    }

    // Temporary storage for subCells
    List<DynamicList<label>> subCells(8);

    scalar deltaT = mesh().time().deltaTValue();

    label collisionCandidates = 0;

    label collisions = 0;

    forAll(cellOccupancy_, celli)
    {
        const DynamicList<ParcelType*>& cellParcels(cellOccupancy_[celli]);

        label nC(cellParcels.size());

        if (nC > 1)
        {
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Assign particles to one of 8 Cartesian subCells

            // Clear temporary lists
            forAll(subCells, i)
            {
                subCells[i].clear();
            }

            // Inverse addressing specifying which subCell a parcel is in
            List<label> whichSubCell(cellParcels.size());

            const point& cC = mesh_.cellCentres()[celli];

            forAll(cellParcels, i)
            {
                const ParcelType& p = *cellParcels[i];
                vector relPos = p.position() - cC;

                label subCell =
                    pos0(relPos.x()) + 2*pos0(relPos.y()) + 4*pos0(relPos.z());

                subCells[subCell].append(i);
                whichSubCell[i] = subCell;
            }

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            scalar sigmaTcRMax = sigmaTcRMax_[celli];

            scalar selectedPairs =
                collisionSelectionRemainder_[celli]
              + 0.5*nC*(nC - 1)*nParticle_*sigmaTcRMax*deltaT
               /mesh_.cellVolumes()[celli];

            label nCandidates(selectedPairs);
            collisionSelectionRemainder_[celli] = selectedPairs - nCandidates;
            collisionCandidates += nCandidates;

            for (label c = 0; c < nCandidates; c++)
            {
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // subCell candidate selection procedure

                // Select the first collision candidate
                label candidateP = rndGen_.position<label>(0, nC - 1);

                // Declare the second collision candidate
                label candidateQ = -1;

                List<label> subCellPs = subCells[whichSubCell[candidateP]];
                label nSC = subCellPs.size();

                if (nSC > 1)
                {
                    // If there are two or more particle in a subCell, choose
                    // another from the same cell.  If the same candidate is
                    // chosen, choose again.

                    do
                    {
                        label i = rndGen_.position<label>(0, nSC - 1);
                        candidateQ = subCellPs[i];
                    } while (candidateP == candidateQ);
                }
                else
                {
                    // Select a possible second collision candidate from the
                    // whole cell.  If the same candidate is chosen, choose
                    // again.

                    do
                    {
                        candidateQ = rndGen_.position<label>(0, nC - 1);
                    } while (candidateP == candidateQ);
                }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                // uniform candidate selection procedure

                // // Select the first collision candidate
                // label candidateP = rndGen_.position<label>(0, nC-1);

                // // Select a possible second collision candidate
                // label candidateQ = rndGen_.position<label>(0, nC-1);

                // // If the same candidate is chosen, choose again
                // while (candidateP == candidateQ)
                // {
                //     candidateQ = rndGen_.position<label>(0, nC-1);
                // }

                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                ParcelType& parcelP = *cellParcels[candidateP];
                ParcelType& parcelQ = *cellParcels[candidateQ];

                scalar sigmaTcR = binaryCollision().sigmaTcR
                (
                    parcelP,
                    parcelQ
                );

                // Update the maximum value of sigmaTcR stored, but use the
                // initial value in the acceptance-rejection criteria because
                // the number of collision candidates selected was based on this

                if (sigmaTcR > sigmaTcRMax_[celli])
                {
                    sigmaTcRMax_[celli] = sigmaTcR;
                }

                if ((sigmaTcR/sigmaTcRMax) > rndGen_.sample01<scalar>())
                {
                    binaryCollision().collide
                    (
                        parcelP,
                        parcelQ
                    );

                    collisions++;
                }
            }
        }
    }

    reduce(collisions, sumOp<label>());

    reduce(collisionCandidates, sumOp<label>());

    sigmaTcRMax_.correctBoundaryConditions();

    if (collisionCandidates)
    {
        Info<< "    Collisions                      = "
            << collisions << nl
            << "    Acceptance rate                 = "
            << scalar(collisions)/scalar(collisionCandidates) << nl
            << endl;
    }
    else
    {
        Info<< "    No collisions" << endl;
    }
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::resetFields()
{
    q_ = dimensionedScalar("0", dimensionSet(1, 0, -3, 0, 0), Zero);
    fD_ = dimensionedVector("0", dimensionSet(1, -1, -2, 0, 0), Zero);

    rhoN_ = dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0), VSMALL);
    rhoM_ =  dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), VSMALL);
    dsmcRhoN_ = dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0), Zero);
    linearKE_ = dimensionedScalar("0", dimensionSet(1, -1, -2, 0, 0), Zero);
    internalE_ = dimensionedScalar("0", dimensionSet(1, -1, -2, 0, 0), Zero);
    iDof_ = dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0), VSMALL);
    momentum_ = dimensionedVector("0", dimensionSet(1, -2, -1, 0, 0), Zero);
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::calculateFields()
{
    scalarField& rhoN = rhoN_.primitiveFieldRef();
    scalarField& rhoM = rhoM_.primitiveFieldRef();
    scalarField& dsmcRhoN = dsmcRhoN_.primitiveFieldRef();
    scalarField& linearKE = linearKE_.primitiveFieldRef();
    scalarField& internalE = internalE_.primitiveFieldRef();
    scalarField& iDof = iDof_.primitiveFieldRef();
    vectorField& momentum = momentum_.primitiveFieldRef();

    for (const ParcelType& p : *this)
    {
        const label celli = p.cell();

        rhoN[celli]++;
        rhoM[celli] += constProps(p.typeId()).mass();
        dsmcRhoN[celli]++;
        linearKE[celli] += 0.5*constProps(p.typeId()).mass()*(p.U() & p.U());
        internalE[celli] += p.Ei();
        iDof[celli] += constProps(p.typeId()).internalDegreesOfFreedom();
        momentum[celli] += constProps(p.typeId()).mass()*p.U();
    }

    rhoN *= nParticle_/mesh().cellVolumes();
    rhoN_.correctBoundaryConditions();

    rhoM *= nParticle_/mesh().cellVolumes();
    rhoM_.correctBoundaryConditions();

    dsmcRhoN_.correctBoundaryConditions();

    linearKE *= nParticle_/mesh().cellVolumes();
    linearKE_.correctBoundaryConditions();

    internalE *= nParticle_/mesh().cellVolumes();
    internalE_.correctBoundaryConditions();

    iDof *= nParticle_/mesh().cellVolumes();
    iDof_.correctBoundaryConditions();

    momentum *= nParticle_/mesh().cellVolumes();
    momentum_.correctBoundaryConditions();
}


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::addNewParcel
(
    const vector& position,
    const label celli,
    const vector& U,
    const scalar Ei,
    const label typeId
)
{
    this->addParticle(new ParcelType(mesh_, position, celli, U, Ei, typeId));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DSMCCloud<ParcelType>::DSMCCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    bool readFields
)
:
    Cloud<ParcelType>(mesh, cloudName, false),
    DSMCBaseCloud(),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(particleProperties_.get<scalar>("nEquivalentParticles")),
    cellOccupancy_(mesh_.nCells()),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    collisionSelectionRemainder_
    (
        IOobject
        (
            IOobject::scopedName(this->name(), "collisionSelectionRemainder"),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    q_
    (
        IOobject
        (
            "q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    fD_
    (
        IOobject
        (
            "fD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    rhoN_
    (
        IOobject
        (
            "rhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    rhoM_
    (
        IOobject
        (
            "rhoM",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    dsmcRhoN_
    (
        IOobject
        (
            "dsmcRhoN",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    linearKE_
    (
        IOobject
        (
            "linearKE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    internalE_
    (
        IOobject
        (
            "internalE",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    iDof_
    (
        IOobject
        (
            "iDof",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    momentum_
    (
        IOobject
        (
            "momentum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    constProps_(),
    rndGen_(Pstream::myProcNo()),
    boundaryT_
    (
        IOobject
        (
            "boundaryT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    boundaryU_
    (
        IOobject
        (
            "boundaryU",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_
    ),
    binaryCollisionModel_
    (
        BinaryCollisionModel<DSMCCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        )
    ),
    wallInteractionModel_
    (
        WallInteractionModel<DSMCCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        )
    ),
    inflowBoundaryModel_
    (
        InflowBoundaryModel<DSMCCloud<ParcelType>>::New
        (
            particleProperties_,
            *this
        )
    )
{
    buildConstProps();
    buildCellOccupancy();

    // Initialise the collision selection remainder to a random value between 0
    // and 1.
    forAll(collisionSelectionRemainder_, i)
    {
        collisionSelectionRemainder_[i] = rndGen_.sample01<scalar>();
    }

    if (readFields)
    {
        ParcelType::readFields(*this);
    }
}


template<class ParcelType>
Foam::DSMCCloud<ParcelType>::DSMCCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const IOdictionary& dsmcInitialiseDict
)
:
    Cloud<ParcelType>(mesh, cloudName, false),
    DSMCBaseCloud(),
    cloudName_(cloudName),
    mesh_(mesh),
    particleProperties_
    (
        IOobject
        (
            cloudName + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    typeIdList_(particleProperties_.lookup("typeIdList")),
    nParticle_(particleProperties_.get<scalar>("nEquivalentParticles")),
    cellOccupancy_(),
    sigmaTcRMax_
    (
        IOobject
        (
            this->name() + "SigmaTcRMax",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, 3, -1, 0, 0), Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    collisionSelectionRemainder_
    (
        IOobject
        (
            IOobject::scopedName(this->name(), "collisionSelectionRemainder"),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    q_
    (
        IOobject
        (
            this->name() + "q_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1, 0, -3, 0, 0), Zero)
    ),
    fD_
    (
        IOobject
        (
            this->name() + "fD_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimensionSet(1, -1, -2, 0, 0), Zero)
    ),
    rhoN_
    (
        IOobject
        (
            this->name() + "rhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0), VSMALL)
    ),
    rhoM_
    (
        IOobject
        (
            this->name() + "rhoM_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(1, -3, 0, 0, 0), VSMALL)
    ),
    dsmcRhoN_
    (
        IOobject
        (
            this->name() + "dsmcRhoN_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(0, -3, 0, 0, 0), Zero)
    ),
    linearKE_
    (
        IOobject
        (
            this->name() + "linearKE_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), Zero)
    ),
    internalE_
    (
        IOobject
        (
            this->name() + "internalE_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), Zero)
    ),
    iDof_
    (
        IOobject
        (
            this->name() + "iDof_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0), VSMALL)
    ),
    momentum_
    (
        IOobject
        (
            this->name() + "momentum_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimensionSet(1, -2, -1, 0, 0), Zero)
    ),
    constProps_(),
    rndGen_(Pstream::myProcNo()),
    boundaryT_
    (
        volScalarField
        (
            IOobject
            (
                "boundaryT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimensionSet(0, 0, 0, 1, 0), Zero)
        )
    ),
    boundaryU_
    (
        volVectorField
        (
            IOobject
            (
                "boundaryU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimensionSet(0, 1, -1, 0, 0), Zero)
        )
    ),
    binaryCollisionModel_(),
    wallInteractionModel_(),
    inflowBoundaryModel_()
{
    clear();
    buildConstProps();
    initialise(dsmcInitialiseDict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DSMCCloud<ParcelType>::~DSMCCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::evolve()
{
    typename ParcelType::trackingData td(*this);

    // Reset the data collection fields
    resetFields();

    if (debug)
    {
        this->dumpParticlePositions();
    }

    // Insert new particles from the inflow boundary
    this->inflowBoundary().inflow();

    // Move the particles ballistically with their current velocities
    Cloud<ParcelType>::move(*this, td, mesh_.time().deltaTValue());

    // Update cell occupancy
    buildCellOccupancy();

    // Calculate new velocities via stochastic collisions
    collisions();

    // Calculate the volume field data
    calculateFields();
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::info() const
{
    label nDSMCParticles = this->size();
    reduce(nDSMCParticles, sumOp<label>());

    scalar nMol = nDSMCParticles*nParticle_;

    vector linearMomentum = linearMomentumOfSystem();
    reduce(linearMomentum, sumOp<vector>());

    scalar linearKineticEnergy = linearKineticEnergyOfSystem();
    reduce(linearKineticEnergy, sumOp<scalar>());

    scalar internalEnergy = internalEnergyOfSystem();
    reduce(internalEnergy, sumOp<scalar>());

    Info<< "Cloud name: " << this->name() << nl
        << "    Number of dsmc particles        = "
        << nDSMCParticles
        << endl;

    if (nDSMCParticles)
    {
        Info<< "    Number of molecules             = "
            << nMol << nl
            << "    Mass in system                  = "
            << returnReduce(massInSystem(), sumOp<scalar>()) << nl
            << "    Average linear momentum         = "
            << linearMomentum/nMol << nl
            << "   |Average linear momentum|        = "
            << mag(linearMomentum)/nMol << nl
            << "    Average linear kinetic energy   = "
            << linearKineticEnergy/nMol << nl
            << "    Average internal energy         = "
            << internalEnergy/nMol << nl
            << "    Average total energy            = "
            << (internalEnergy + linearKineticEnergy)/nMol
            << endl;
    }
}


template<class ParcelType>
Foam::vector Foam::DSMCCloud<ParcelType>::equipartitionLinearVelocity
(
    scalar temperature,
    scalar mass
)
{
    return
        sqrt(physicoChemical::k.value()*temperature/mass)
       *rndGen_.GaussNormal<vector>();
}


template<class ParcelType>
Foam::scalar Foam::DSMCCloud<ParcelType>::equipartitionInternalEnergy
(
    scalar temperature,
    direction iDof
)
{
    if (iDof == 0)
    {
        return 0;
    }
    else if (iDof == 2)
    {
        // Special case for iDof = 2, i.e. diatomics;
        return
        (
            -log(rndGen_.sample01<scalar>())
            *physicoChemical::k.value()*temperature
        );
    }


    const scalar a = 0.5*iDof - 1;
    scalar energyRatio = 0;
    scalar P = -1;

    do
    {
        energyRatio = 10*rndGen_.sample01<scalar>();
        P = pow((energyRatio/a), a)*exp(a - energyRatio);
    } while (P < rndGen_.sample01<scalar>());

    return energyRatio*physicoChemical::k.value()*temperature;
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::dumpParticlePositions() const
{
    OFstream pObj
    (
        this->db().time().path()/"parcelPositions_"
      + this->name() + "_"
      + this->db().time().timeName() + ".obj"
    );

    for (const ParcelType& p : *this)
    {
        pObj<< "v " << p.position().x()
            << ' '  << p.position().y()
            << ' '  << p.position().z()
            << nl;
    }

    pObj.flush();
}


template<class ParcelType>
void Foam::DSMCCloud<ParcelType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<ParcelType>::autoMap(mapper);

    // Update the cell occupancy field
    cellOccupancy_.setSize(mesh_.nCells());
    buildCellOccupancy();

    // Update the inflow BCs
    this->inflowBoundary().autoMap(mapper);
}


// ************************************************************************* //

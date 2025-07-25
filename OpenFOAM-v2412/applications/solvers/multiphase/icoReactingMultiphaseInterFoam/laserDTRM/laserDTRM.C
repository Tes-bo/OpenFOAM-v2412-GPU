/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "laserDTRM.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "constants.H"
#include "unitConversion.H"
#include "interpolationCell.H"
#include "interpolationCellPoint.H"
#include "Random.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(laserDTRM, 0);
        addToRadiationRunTimeSelectionTables(laserDTRM);
    }

    defineTemplateTypeNameAndDebugWithName
    (
        Cloud<DTRMParticle>,
        "DTRMCloud",
        0
    );

} // End namespace Foam


const Foam::Enum<Foam::radiation::laserDTRM::powerDistributionMode>
Foam::radiation::laserDTRM::powerDistNames_
{
    { powerDistributionMode::pdGaussian, "Gaussian" },
    { powerDistributionMode::pdManual,   "manual" },
    { powerDistributionMode::pdUniform,  "uniform" },
    { powerDistributionMode::pdGaussianPeak, "GaussianPeak" },
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::radiation::laserDTRM::calculateIp(scalar r, scalar theta)
{
    const scalar t = mesh_.time().value();
    const scalar power = laserPower_->value(t);
    switch (mode_)
    {
        case pdGaussianPeak:
        {
            return I0_*exp(-2.0*sqr(r)/sqr(sigma_));
            break;
        }
        case pdGaussian:
        {
            scalar I0 = power/(mathematical::twoPi*sqr(sigma_));

            return I0*exp(-sqr(r)/2.0/sqr(sigma_));
            break;
        }
        case pdManual:
        {
            return power*powerDistribution_()(theta, r);
            break;
        }
        case pdUniform:
        {
            return power/(mathematical::pi*sqr(focalLaserRadius_));
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled type " << powerDistNames_[mode_]
                << abort(FatalError);
            break;
        }
    }

    return 0;
}


Foam::tmp<Foam::volVectorField> Foam::radiation::laserDTRM::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-7/cbrt(average(mesh_.V()))
    );

    const volVectorField gradAlphaf
    (
         alpha2*fvc::grad(alpha1)
       - alpha1*fvc::grad(alpha2)
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf)+ deltaN);
}


void Foam::radiation::laserDTRM::initialiseReflection()
{
    if (found("reflectionModel"))
    {
        dictTable modelDicts(lookup("reflectionModel"));

        forAllConstIters(modelDicts, iter)
        {
            const phasePairKey& key = iter.key();

            reflections_.insert
            (
                key,
                reflectionModel::New
                (
                    iter.val(),
                    mesh_
                )
            );
        }

        reflectionSwitch_ = returnReduceOr(reflections_.size());
    }
}


void Foam::radiation::laserDTRM::initialise()
{
    // Initialise the DTRM particles
    DTRMCloud_.clear();

    const scalar t = mesh_.time().value();
    const vector lPosition = focalLaserPosition_->value(t);
    const vector lDir = normalised(laserDirection_->value(t));

    DebugInfo
        << "Laser position : " << lPosition << nl
        << "Laser direction : " << lDir << endl;

    // Find a vector on the area plane. Normal to laser direction
    vector rArea = Zero;
    scalar magr = 0.0;

    {
        Random rnd(1234);

        while (magr < VSMALL)
        {
            vector v = rnd.sample01<vector>();
            rArea = v - (v & lDir)*lDir;
            magr = mag(rArea);
        }
    }
    rArea.normalise();

    scalar dr =  focalLaserRadius_/ndr_;
    scalar dTheta =  mathematical::twoPi/ndTheta_;

    nParticles_ = ndr_*ndTheta_;

    switch (mode_)
    {
        case pdGaussianPeak:
        {
            I0_ = get<scalar>("I0");
            sigma_ = get<scalar>("sigma");
            break;
        }
        case pdGaussian:
        {
            sigma_ = get<scalar>("sigma");
            break;
        }
        case pdManual:
        {
            powerDistribution_.reset
            (
                new interpolation2DTable<scalar>(*this)
            );

            break;
        }
        case pdUniform:
        {
            break;
        }
    }

    // Count the number of missed positions
    label nMissed = 0;

    // Target position
    point p1 = vector::zero;

    // Seed DTRM particles
    // TODO: currently only applicable to 3-D cases
    point p0 = lPosition;
    scalar power(0.0);
    scalar area(0.0);
    p1 = p0;
    if (mesh_.nGeometricD() == 3)
    {
        for (label ri = 0; ri < ndr_; ri++)
        {
            scalar r1 = SMALL + dr*ri;

            scalar r2 = r1 + dr;

            scalar rP = ((r1 + r2)/2);

            // local radius on disk
            vector localR = ((r1 + r2)/2)*rArea;

            // local final radius on disk
            vector finalR = rP*rArea;

            scalar theta0 = 0.0;//dTheta/2.0;
            for (label thetai = 0; thetai < ndTheta_; thetai++)
            {
                scalar theta1 = theta0 + SMALL  + dTheta*thetai;

                scalar theta2 = theta1 + dTheta;

                scalar thetaP = (theta1 + theta2)/2.0;

                quaternion Q(lDir, thetaP);

                // Initial position on disk
                vector initialPos = (Q.R() & localR);

                // Final position on disk
                vector finalPos = (Q.R() & finalR);

                // Initial position
                p0 = lPosition + initialPos;

                // calculate target point using new deviation rl
                p1 = lPosition + finalPos + (0.5*maxTrackLength_*lDir);

                scalar Ip = calculateIp(rP, thetaP);

                scalar dAi = (sqr(r2) - sqr(r1))*(theta2 - theta1)/2.0;

                power += Ip*dAi;
                area += dAi;

                label cellI = mesh_.findCell(p0);

                if (cellI != -1)
                {
                    // Create a new particle
                    DTRMParticle* pPtr =
                        new DTRMParticle(mesh_, p0, p1, Ip, cellI, dAi, -1);

                    // Add to cloud
                    DTRMCloud_.addParticle(pPtr);
                }

                if (nMissed < 10 && returnReduceAnd(cellI < 0))
                {
                    ++nMissed;
                    WarningInFunction
                        << "Cannot find owner cell for focalPoint at "
                        << p0 << endl;
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Current functionality limited to 3-D cases"
            << exit(FatalError);
    }

    if (nMissed)
    {
        Info<< "Seeding missed " << nMissed << " locations" << endl;
    }

    DebugInfo
        << "Total Power in the laser : " << power << nl
        << "Total Area in the laser : " << area << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::laserDTRM::laserDTRM(const volScalarField& T)
:
    radiationModel(typeName, T),
    mode_(powerDistNames_.get("mode", *this)),
    DTRMCloud_(mesh_, Foam::zero{}, "DTRMCloud"),  // Empty cloud
    nParticles_(0),
    ndTheta_(get<label>("nTheta")),
    ndr_(get<label>("nr")),
    maxTrackLength_(mesh_.bounds().mag()),

    focalLaserPosition_
    (
        Function1<point>::New("focalLaserPosition", *this, &mesh_)
    ),

    laserDirection_
    (
        Function1<vector>::New("laserDirection", *this, &mesh_)
    ),

    focalLaserRadius_(get<scalar>("focalLaserRadius")),
    qualityBeamLaser_
    (
        getOrDefault<scalar>("qualityBeamLaser", 0)
    ),

    sigma_(0),
    I0_(0),
    laserPower_(Function1<scalar>::New("laserPower", *this, &mesh_)),
    powerDistribution_(),

    reflectionSwitch_(false),

    alphaCut_(getOrDefault<scalar>("alphaCut", 0.5)),

    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    ),
    Q_
    (
        IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimPower/dimVolume, Zero)
    )
{
    initialiseReflection();

    initialise();
}


Foam::radiation::laserDTRM::laserDTRM
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    mode_(powerDistNames_.get("mode", *this)),
    DTRMCloud_(mesh_, Foam::zero{}, "DTRMCloud"),  // Empty cloud
    nParticles_(0),
    ndTheta_(get<label>("nTheta")),
    ndr_(get<label>("nr")),
    maxTrackLength_(mesh_.bounds().mag()),

    focalLaserPosition_
    (
        Function1<point>::New("focalLaserPosition", *this, &mesh_)
    ),
    laserDirection_
    (
        Function1<vector>::New("laserDirection", *this, &mesh_)
    ),

    focalLaserRadius_(get<scalar>("focalLaserRadius")),
    qualityBeamLaser_
    (
        getOrDefault<scalar>("qualityBeamLaser", 0)
    ),

    sigma_(0),
    I0_(0),
    laserPower_(Function1<scalar>::New("laserPower", *this, &mesh_)),
    powerDistribution_(),

    reflectionSwitch_(false),

    alphaCut_(getOrDefault<scalar>("alphaCut", 0.5)),

    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    e_
    (
        IOobject
        (
            "e",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimLength/pow3(dimTime), Zero)
    ),
    Q_
    (
        IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimPower/pow3(dimLength), Zero)
    )
{
    initialiseReflection();
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::laserDTRM::read()
{
    if (radiationModel::read())
    {
        return true;
    }

    return false;
}

Foam::label Foam::radiation::laserDTRM::nBands() const
{
    return 1;
}


void Foam::radiation::laserDTRM::calculate()
{
    auto treflectingCells = volScalarField::New
    (
        "reflectingCellsVol",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar("zero", dimless, -1)
    );
    auto& reflectingCellsVol = treflectingCells.ref();

    auto tnHat = volVectorField::New
    (
        "nHat",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedVector(dimless, Zero)
    );
    auto& nHat = tnHat.ref();


    // Reset the field
    Q_ == dimensionedScalar(Q_.dimensions(), Zero);

    a_ = absorptionEmission_->a();
    e_ = absorptionEmission_->e();
    E_ = absorptionEmission_->E();

    const interpolationCell<scalar> aInterp(a_);
    const interpolationCell<scalar> eInterp(e_);
    const interpolationCell<scalar> EInterp(E_);
    const interpolationCell<scalar> TInterp(T_);

    labelField reflectingCells(mesh_.nCells(), -1);

    UPtrList<reflectionModel> reflectionUPtr;

    if (reflectionSwitch_)
    {
        reflectionUPtr.resize(reflections_.size());

        label reflectionModelId(0);
        forAllIters(reflections_, iter1)
        {
            reflectionModel& model = iter1()();

            reflectionUPtr.set(reflectionModelId, &model);

            const volScalarField& alphaFrom =
                mesh_.lookupObject<volScalarField>
                (
                    IOobject::groupName("alpha", iter1.key().first())
                );

            const volScalarField& alphaTo =
                mesh_.lookupObject<volScalarField>
                (
                    IOobject::groupName("alpha", iter1.key().second())
                );

            const volVectorField nHatPhase(nHatfv(alphaFrom, alphaTo));

            const volScalarField gradAlphaf
            (
                fvc::grad(alphaFrom)
              & fvc::grad(alphaTo)
            );

            const volScalarField nearInterface(pos(alphaTo - alphaCut_));

            const volScalarField mask(nearInterface*gradAlphaf);

            forAll(alphaFrom, cellI)
            {
                if
                (
                    nearInterface[cellI]
                 && mag(nHatPhase[cellI]) > 0.99
                 && mask[cellI] < 0
                )
                {
                    reflectingCells[cellI] = reflectionModelId;
                    reflectingCellsVol[cellI] = reflectionModelId;
                    if (mag(nHat[cellI]) == 0.0)
                    {
                        nHat[cellI] += nHatPhase[cellI];
                    }
                }
            }
            reflectionModelId++;
        }
    }

    interpolationCellPoint<vector> nHatInterp(nHat);

    DTRMParticle::trackingData td
    (
        DTRMCloud_,
        aInterp,
        eInterp,
        EInterp,
        TInterp,
        nHatInterp,
        reflectingCells,
        reflectionUPtr,
        Q_
    );

    Info<< "Move particles..."
        << returnReduce(DTRMCloud_.size(), sumOp<label>()) << endl;

    DTRMCloud_.move(DTRMCloud_, td, mesh_.time().deltaTValue());

    // Normalize by cell volume
    Q_.primitiveFieldRef() /= mesh_.V();

    if (debug)
    {
        Info<< "Final number of particles..."
            << returnReduce(DTRMCloud_.size(), sumOp<label>()) << endl;

        pointField lines(2*DTRMCloud_.size());
        {
            label i = 0;
            for (const DTRMParticle& p : DTRMCloud_)
            {
                lines[i] = p.position();
                lines[i+1] = p.p0();
                i += 2;
            }
        }

        globalIndex::gatherInplaceOp(lines);

        if (UPstream::master())
        {
            OBJstream os(type() + "-particlePath.obj");

            for (label pointi = 0; pointi < lines.size(); pointi += 2)
            {
                os.writeLine(lines[pointi], lines[pointi+1]);
            }
        }

        scalar totalQ = gSum(Q_.primitiveFieldRef()*mesh_.V());
        Info << "Total energy absorbed [W]: " << totalQ << endl;

        if (mesh_.time().writeTime())
        {
            reflectingCellsVol.write();
            nHat.write();
        }
    }

    // Clear and initialise the cloud
    // NOTE: Possible to reset original particles, but this requires
    // data transfer for the cloud in differet processors.
    initialise();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::laserDTRM::Rp() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "zero",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimPower/dimVolume/pow4(dimTemperature), Zero)
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::laserDTRM::Ru() const
{
    return Q_.internalField();
}


// ************************************************************************* //

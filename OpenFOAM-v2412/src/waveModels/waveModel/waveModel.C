/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "waveModel.H"
#include "fvMesh.H"
#include "polyPatch.H"
#include "gravityMeshObject.H"
#include "volFields.H"
#include "fvPatchFields.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(waveModel, 0);
    defineRunTimeSelectionTable(waveModel, patch);
}

const Foam::word Foam::waveModel::dictName("waveProperties");


Foam::word Foam::waveModel::modelName(const word& patchName)
{
    return dictName + '.' + patchName;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::waveModel::initialiseGeometry()
{
    // Determine local patch coordinate system given by:
    // - X: streamwise: patch normal
    // - Y: spanwise: Z^X
    // - Z: up: (negative) gravity direction
    vector x = normalised(-gAverage(patch_.faceAreas()));
    vector z = -g_/mag(g_);
    vector y = z^x;

    // Rotation from global<->local about global origin
    Rlg_ = tensor(x, y, z);
    Rgl_ = Rlg_.T();

    // Global patch extents
    const vectorField& Cp = patch_.localPoints();
    const vectorField CpLocal(Rgl_ & Cp);
    boundBox bb(CpLocal, true);
    const scalar xMin = bb.min().x();
    const scalar xMax = bb.max().x();
    const scalar yMin = bb.min().y();
    const scalar yMax = bb.max().y();
    zSpan_ = bb.max().z() - bb.min().z();

    // Global x, y positions of the paddle centres
    xPaddle_.setSize(nPaddle_, 0);
    yPaddle_.setSize(nPaddle_, 0);
    const scalar xMid = xMin + 0.5*(xMax - xMin);
    const scalar paddleDy = (yMax - yMin)/scalar(nPaddle_);
    for (label paddlei = 0; paddlei < nPaddle_; ++paddlei)
    {
        xPaddle_[paddlei] = xMid;
        yPaddle_[paddlei] = paddlei*paddleDy + yMin + 0.5*paddleDy;
    }

    // Local face centres
    const vectorField CfLocal(Rgl_ & patch_.faceCentres());
    z_ = CfLocal.component(2);

    // Local face extents in z-direction
    zMin_.setSize(patch_.size());
    zMax_.setSize(patch_.size());
    const faceList& faces = patch_.localFaces();
    forAll(faces, facei)
    {
        const face& f = faces[facei];
        const label nPoint = f.size();
        zMin_[facei] = CpLocal[f[0]].z();
        zMax_[facei] = CpLocal[f[0]].z();

        for (label fpi = 1; fpi < nPoint; ++fpi)
        {
            const label pointi = f[fpi];
            zMin_[facei] = min(zMin_[facei], CpLocal[pointi].z());
            zMax_[facei] = max(zMax_[facei], CpLocal[pointi].z());
        }
    }

    // Set minimum z reference level
    zMin0_ = gMin(zMin_);

    // Local paddle-to-face addressing
    faceToPaddle_.setSize(patch_.size(), -1);
    forAll(faceToPaddle_, facei)
    {
        faceToPaddle_[facei] = floor((CfLocal[facei].y() - yMin)/paddleDy);
    }
}


Foam::tmp<Foam::scalarField> Foam::waveModel::waterLevel() const
{
    // Note: initialising as initial depth
    auto tlevel = tmp<scalarField>::New(nPaddle_, initialDepth_);
    auto& level = tlevel.ref();

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);
    const fvPatchScalarField& alphap = alpha.boundaryField()[patch_.index()];
    const scalarField alphac(alphap.patchInternalField());

    const scalarField& magSf = alphap.patch().magSf();
    scalarList paddleMagSf(nPaddle_, Zero);
    scalarList paddleWettedMagSf(nPaddle_, Zero);

    forAll(alphac, facei)
    {
        label paddlei = faceToPaddle_[facei];
        paddleMagSf[paddlei] += magSf[facei];
        paddleWettedMagSf[paddlei] += magSf[facei]*alphac[facei];
    }

    forAll(paddleMagSf, paddlei)
    {
        reduce(paddleMagSf[paddlei], sumOp<scalar>());
        reduce(paddleWettedMagSf[paddlei], sumOp<scalar>());
        level[paddlei] +=
            paddleWettedMagSf[paddlei]*zSpan_
           /(paddleMagSf[paddlei] + ROOTVSMALL);
    }

    return tlevel;
}


void Foam::waveModel::setAlpha(const scalarField& level)
{
    forAll(alpha_, facei)
    {
        const label paddlei = faceToPaddle_[facei];
        const scalar paddleCalc = level[paddlei];

        const scalar zMin0 = zMin_[facei] - zMin0_;
        const scalar zMax0 = zMax_[facei] - zMin0_;

        if (zMax0 < paddleCalc)
        {
            alpha_[facei] = 1.0;
        }
        else if (zMin0 > paddleCalc)
        {
            alpha_[facei] = 0.0;
        }
        else
        {
            scalar dz = paddleCalc - zMin0;
            alpha_[facei] = dz/(zMax0 - zMin0);
        }
    }
}


void Foam::waveModel::setPaddlePropeties
(
    const scalarField& level,
    const label facei,
    scalar& fraction,
    scalar& z
) const
{
    const label paddlei = faceToPaddle_[facei];
    const scalar paddleCalc = level[paddlei];
    const scalar paddleHeight = min(paddleCalc, waterDepthRef_);
    const scalar zMin = zMin_[facei] - zMin0_;
    const scalar zMax = zMax_[facei] - zMin0_;

    fraction = 1;
    z = 0;

    if (zMax < paddleHeight)
    {
        z = z_[facei] - zMin0_;
    }
    else if (zMin > paddleCalc)
    {
        fraction = -1;
    }
    else
    {
        if (paddleCalc < waterDepthRef_)
        {
            if ((zMax > paddleCalc) && (zMin < paddleCalc))
            {
                scalar dz = paddleCalc - zMin;
                fraction = dz/(zMax - zMin);
                z = z_[facei] - zMin0_;
            }
        }
        else
        {
            if (zMax < paddleCalc)
            {
                z = waterDepthRef_;
            }
            else if ((zMax > paddleCalc) && (zMin < paddleCalc))
            {
                scalar dz = paddleCalc - zMin;
                fraction = dz/(zMax - zMin);
                z = waterDepthRef_;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModel::waveModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    IOdictionary
    (
        IOobject
        (
            modelName(patch.name()),
            Time::timeName(mesh.time().startTime().value()),
            "uniform",
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    mesh_(mesh),
    patch_(patch),
    g_(meshObjects::gravity::New(mesh.time()).value()),
    UName_("U"),
    alphaName_("alpha"),
    Rgl_(tensor::I),
    Rlg_(tensor::I),
    nPaddle_(1),
    xPaddle_(),
    yPaddle_(),
    z_(),
    zSpan_(0),
    zMin_(),
    zMax_(),
    waterDepthRef_(0),
    initialDepth_(0),
    currTimeIndex_(-1),
    activeAbsorption_(false),
    U_(patch.size(), Zero),
    alpha_(patch.size(), Zero)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModel::readDict(const dictionary& overrideDict)
{
    readOpt(IOobject::READ_IF_PRESENT);
    if (headerOk())
    {
        IOdictionary::regIOobject::read();
    }

    merge(overrideDict);

    readIfPresent("U", UName_);
    readIfPresent("alpha", alphaName_);

    readEntry("nPaddle", nPaddle_);
    if (nPaddle_ < 1)
    {
        FatalIOErrorInFunction(*this)
            << "Number of paddles must be greater than zero.  Supplied"
            << " value nPaddles = " << nPaddle_
            << exit(FatalIOError);
    }

    readIfPresent("initialDepth", initialDepth_);

    // Need to initialise the geometry before calling waterLevel()
    initialiseGeometry();

    // Set the reference water depth
    if (!readIfPresent("waterDepthRef", waterDepthRef_))
    {
        scalar waterDepth = 0;
        if (readIfPresent("waterDepth", waterDepth))
        {
            waterDepthRef_ = waterDepth;
        }
        else
        {
            const scalarField level(waterLevel());
            if (level.size())
            {
                waterDepthRef_ = level.first();
            }
        }

        // Avoid potential zero...
        waterDepthRef_ += SMALL;

        // Insert the reference water depth into [this] to enable restart
        add("waterDepthRef", waterDepthRef_);
    }

    return true;
}


void Foam::waveModel::correct(const scalar t)
{
    if (mesh_.time().timeIndex() != currTimeIndex_)
    {
        Info<< "Updating " << type() << " wave model for patch "
            << patch_.name() << endl;

        // Time ramp weight
        const scalar tCoeff = timeCoeff(t);

        // Reset the velocity and phase fraction fields
        U_ = vector::zero;
        alpha_ = 0;

        // Update the calculated water level field
        scalarField calculatedLevel(nPaddle_, Zero);

        if (patch_.size())
        {
            // Set wave level
            setLevel(t, tCoeff, calculatedLevel);

            // Update the velocity field
            setVelocity(t, tCoeff, calculatedLevel);

            // Update the phase fraction field
            setAlpha(calculatedLevel);
        }

        if (activeAbsorption_)
        {
            const scalarField activeLevel(this->waterLevel());

            forAll(U_, facei)
            {
                const label paddlei = faceToPaddle_[facei];

                if (zMin_[facei] - zMin0_ < activeLevel[paddlei])
                {
                    scalar UCorr =
                        (calculatedLevel[paddlei] - activeLevel[paddlei])
                       *sqrt(mag(g_)/activeLevel[paddlei]);

                    U_[facei].x() += UCorr;
                }
                else
                {
                    U_[facei].x() = 0;
                }
            }
        }

        // Transform velocity into global coordinate system
        U_ = Rlg_ & U_;

        currTimeIndex_ = mesh_.time().timeIndex();
    }
}


const Foam::vectorField& Foam::waveModel::U() const
{
    return U_;
}


const Foam::scalarField& Foam::waveModel::alpha() const
{
    return alpha_;
}


void Foam::waveModel::info(Ostream& os) const
{
    os  << "Wave model: patch " << patch_.name() << nl
        << "    Type : " << type() << nl
        << "    Velocity field name : " << UName_ << nl
        << "    Phase fraction field name : " << alphaName_ << nl
        << "    Transformation from local to global system : " << Rlg_ << nl
        << "    Number of paddles: " << nPaddle_ << nl
        << "    Reference water depth : " << waterDepthRef_ << nl
        << "    Active absorption: " << activeAbsorption_ << nl;
}


// ************************************************************************* //

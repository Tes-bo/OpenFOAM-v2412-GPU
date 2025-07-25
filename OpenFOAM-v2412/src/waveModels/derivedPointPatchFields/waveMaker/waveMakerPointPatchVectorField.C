/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 IH-Cantabria
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "waveMakerPointPatchVectorField.H"
#include "mathematicalConstants.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "gravityMeshObject.H"

#include "polyMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum<Foam::waveMakerPointPatchVectorField::motionTypes>
Foam::waveMakerPointPatchVectorField::motionTypeNames
({
    { motionTypes::piston, "piston" },
    { motionTypes::flap, "flap" },
    { motionTypes::solitary, "solitary" }
});


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const Foam::vector& Foam::waveMakerPointPatchVectorField::g()
{
    const meshObjects::gravity& gf = meshObjects::gravity::New(db().time());

    if (mag(gf.value()) < SMALL)
    {
        FatalErrorInFunction
            << "Gravity vector is not set.  Please update "
            << gf.uniformDimensionedVectorField::path()
            << exit(FatalError);
    }

    return gf.value();
}


Foam::scalar Foam::waveMakerPointPatchVectorField::waveLength
(
    const scalar h,
    const scalar T
)
{
    const scalar L0 = mag(g())*T*T/(constant::mathematical::twoPi);
    scalar L = L0;

    for (label i=1; i<=100; ++i)
    {
        L = L0*tanh(constant::mathematical::twoPi*h/L);
    }

    return L;
}


Foam::scalar Foam::waveMakerPointPatchVectorField::timeCoeff
(
    const scalar t
) const
{
    return clamp(t/rampTime_, zero_one{});
}


void Foam::waveMakerPointPatchVectorField::initialiseGeometry()
{
    // Global patch extents
    const vectorField& Cp = this->patch().localPoints();
    const vectorField CpLocal(Cp);
    boundBox bb(CpLocal, true);

    const scalar xMin = bb.min().x();
    const scalar xMax = bb.max().x();
    const scalar yMin = bb.min().y();
    const scalar yMax = bb.max().y();
    zSpan_ = bb.max().z() - bb.min().z();

    zMinGb_ = bb.min().z();
    reduce(zMinGb_, minOp<scalar>());

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
    x_ = this->patch().localPoints().component(0);
    y_ = this->patch().localPoints().component(1);
    z_ = this->patch().localPoints().component(2);

    // Local point-to-paddle addressing
    pointToPaddle_.setSize(this->patch().size(), -1);

    forAll(pointToPaddle_, ppi)
    {
        pointToPaddle_[ppi] = floor((y_[ppi] - yMin)/(paddleDy+0.01*paddleDy));
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    motionType_(motionTypes::piston),
    n_(Zero),
    gHat_(Zero),
    initialDepth_(0),
    wavePeriod_(0),
    waveHeight_(0),
    wavePhase_(0),
    waveAngle_(0),
    startTime_(0),
    rampTime_(1),
    secondOrder_(false),
    nPaddle_(0)
{}


Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict, IOobjectOption::NO_READ),
    motionType_(motionTypeNames.get("motionType", dict)),
    n_(dict.get<vector>("n")),
    gHat_(Zero),
    initialDepth_(dict.get<scalar>("initialDepth")),
    wavePeriod_(dict.get<scalar>("wavePeriod")),
    waveHeight_(dict.get<scalar>("waveHeight")),
    wavePhase_(dict.get<scalar>("wavePhase")),
    waveAngle_(dict.getOrDefault<scalar>("waveAngle", 0)),
    startTime_
    (
        dict.getOrDefault<scalar>
        (
            "startTime",
            db().time().startTime().value()
        )
    ),
    rampTime_(dict.get<scalar>("rampTime")),
    secondOrder_(dict.getOrDefault<bool>("secondOrder", false)),
    nPaddle_(dict.getOrDefault<label>("nPaddle", 1))
{
    // Create the co-ordinate system
    if (mag(n_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Patch normal direction vector is not set. 'n' = " << n_
            << exit(FatalIOError);
    }
    n_.normalise();

    gHat_ = (g() - n_*(n_&g()));
    if (mag(gHat_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Patch normal and gravity directions must not be aligned. "
            << "'n' = " << n_ << " 'g' = " << g()
            << exit(FatalIOError);
    }
    gHat_.normalise();

    waveAngle_ *= constant::mathematical::pi/180;

    initialiseGeometry();

    waterDepthRef_.setSize(nPaddle_, -1);

    if (!dict.found("value"))
    {
        updateCoeffs();
    }
}


Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const waveMakerPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    motionType_(ptf.motionType_),
    n_(ptf.n_),
    gHat_(ptf.gHat_),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    wavePhase_(ptf.wavePhase_),
    waveAngle_(ptf.waveAngle_),
    startTime_(ptf.startTime_),
    rampTime_(ptf.rampTime_),
    secondOrder_(ptf.secondOrder_),
    nPaddle_(ptf.nPaddle_)
{}


Foam::waveMakerPointPatchVectorField::waveMakerPointPatchVectorField
(
    const waveMakerPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    motionType_(ptf.motionType_),
    n_(ptf.n_),
    gHat_(ptf.gHat_),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    wavePhase_(ptf.wavePhase_),
    waveAngle_(ptf.waveAngle_),
    startTime_(ptf.startTime_),
    rampTime_(ptf.rampTime_),
    secondOrder_(ptf.secondOrder_),
    nPaddle_(ptf.nPaddle_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveMakerPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (firstTime == 0)
    {
        // Set the reference water depth
        if (initialDepth_ != 0 )
        {
            forAll(waterDepthRef_, paddlei)
            {
                waterDepthRef_[paddlei] = initialDepth_;
            }
        }
        else
        {
            FatalErrorInFunction
               << "initialDepth is not set.  Please update "
               << abort(FatalError);
        }


        Info<< " WaterDepth at the wavepaddles = " << waterDepthRef_ << endl;
        firstTime = 1;
    }

    const scalar t = db().time().value() - startTime_;

    scalarField waveLength_(nPaddle_, -1);

    scalarField waveK(nPaddle_, -1);
    scalarField waveKx(nPaddle_, -1);
    scalarField waveKy(nPaddle_, -1);

    forAll(waveK, padddlei)
    {
        waveLength_[padddlei] =
            waveLength(waterDepthRef_[padddlei], wavePeriod_);

        waveK[padddlei] = constant::mathematical::twoPi/waveLength_[padddlei];
        waveKx[padddlei] = waveK[padddlei]*cos(waveAngle_);
        waveKy[padddlei] = waveK[padddlei]*sin(waveAngle_);
    }
    const scalar sigma = 2*constant::mathematical::pi/wavePeriod_;

    switch (motionType_)
    {
        case motionTypes::flap:
        {
            const pointField& points = patch().localPoints();
            scalarField motionX(patch().localPoints().size(), -1);

            forAll(points, pointi)
            {
                const label paddlei = pointToPaddle_[pointi];

                const scalar phaseTot =
                   waveKx[paddlei]*xPaddle_[paddlei]
                 + waveKy[paddlei]*yPaddle_[paddlei];

                const scalar depthRef = waterDepthRef_[paddlei];
                const scalar kh = waveK[paddlei]*depthRef;
                const scalar pz = points[pointi].component(2);

                const scalar m1 =
                    (4*sinh(kh)/(sinh(2*kh) + 2*kh))
                  * (sinh(kh) + 1/kh*(1 - cosh(kh)));

                const scalar boardStroke = waveHeight_/m1;

                motionX[pointi] = 0.5*boardStroke*sin(phaseTot - sigma*t);

                if (secondOrder_)
                {
                    motionX[pointi] +=
                        sqr(waveHeight_)/(16*depthRef)
                      * (3*cosh(kh)/pow3(sinh(kh)) - 2/m1)
                      * sin(phaseTot - 2*sigma*t);

                }

                motionX[pointi] *= 1.0 + (pz - zMinGb_ - depthRef)/depthRef;

            }

            Field<vector>::operator=(timeCoeff(t)*n_*motionX);

            break;
        }
        case motionTypes::piston:
        {
            const pointField& points = patch().localPoints();
            scalarField motionX(patch().localPoints().size(), -1);

            forAll(points, pointi)
            {
                const label paddlei = pointToPaddle_[pointi];

                const scalar phaseTot =
                    waveKx[paddlei]*xPaddle_[paddlei]
                  + waveKy[paddlei]*yPaddle_[paddlei];

                const scalar depthRef = waterDepthRef_[paddlei];
                const scalar kh = waveK[paddlei]*depthRef;
                const scalar m1 = 2*(cosh(2*kh) - 1.0)/(sinh(2*kh) + 2*kh);
                const scalar boardStroke = waveHeight_/m1;

                motionX[pointi] = 0.5*boardStroke*sin(phaseTot - sigma*t);

                if (secondOrder_)
                {
                    motionX[pointi] +=
                      + sqr(waveHeight_)
                      / (32*depthRef)*(3*cosh(kh)/pow3(sinh(kh)) - 2.0/m1)
                      * sin(phaseTot - 2*sigma*t);
                }
            }

            Field<vector>::operator=(timeCoeff(t)*n_*motionX);

            break;
        }
        case motionTypes::solitary:
        {
            const pointField& points = patch().localPoints();
            scalarField motionX(patch().localPoints().size(), -1);
            const scalar magG = mag(g());

            forAll(points, pointi)
            {
                const label paddlei = pointToPaddle_[pointi];
                const scalar depthRef = waterDepthRef_[paddlei];

                const scalar kappa = sqrt(0.75*waveHeight_/pow3(depthRef));
                const scalar celerity = sqrt(magG*(depthRef + waveHeight_));
                const scalar stroke = sqrt(16*waveHeight_*depthRef/3.0);
                const scalar hr = waveHeight_/depthRef;
                wavePeriod_ = 2.0/(kappa*celerity)*(3.8 + hr);
                const scalar tSolitary = -0.5*wavePeriod_ + t;

                // Newton-Raphson
                scalar theta1 = 0;
                scalar theta2 = 0;
                scalar er = 10000;
                const scalar error = 0.001;
                while (er > error)
                {
                    theta2 =
                        theta1
                      - (theta1 - kappa*celerity*tSolitary + hr*tanh(theta1))
                       /(1.0 + hr*(1.0/cosh(theta1))*(1.0/cosh(theta1)));

                    er = mag(theta1 - theta2);
                    theta1 = theta2;
                }

                motionX[pointi] =
                    waveHeight_/(kappa*depthRef)*tanh(theta1) + 0.5*stroke;
            }

            Field<vector>::operator=(n_*motionX);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << motionTypeNames[motionType_]
                << abort(FatalError);
        }
    }

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::waveMakerPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("motionType", motionTypeNames[motionType_]);
    os.writeEntry("n", n_);
    os.writeEntry("initialDepth", initialDepth_);
    os.writeEntry("wavePeriod", wavePeriod_);
    os.writeEntry("waveHeight", waveHeight_);
    os.writeEntry("wavePhase", wavePhase_);
    os.writeEntry("waveAngle", waveAngle_);
    os.writeEntry("startTime", startTime_);
    os.writeEntry("rampTime", rampTime_);
    os.writeEntry("secondOrder", secondOrder_);
    os.writeEntry("nPaddle", nPaddle_);
    this->writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        waveMakerPointPatchVectorField
    );
}

// ************************************************************************* //

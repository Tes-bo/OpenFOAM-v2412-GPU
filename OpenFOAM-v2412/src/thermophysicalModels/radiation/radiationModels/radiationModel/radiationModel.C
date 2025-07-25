/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"
#include "fvmSup.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(radiationModel, 0);
        defineRunTimeSelectionTable(radiationModel, T);
        defineRunTimeSelectionTable(radiationModel, dictionary);
    }
}

const Foam::word Foam::radiation::radiationModel::externalRadHeatFieldName_ =
    "qrExt";

const Foam::word Foam::radiation::radiationModel::primaryFluxName_ =
    "qprimaryRad";

 const Foam::word Foam::radiation::radiationModel::relfectedFluxName_ =
    "qreflective";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiation::radiationModel::createIOobject
(
    const fvMesh& mesh
) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
    }
    else
    {
        io.readOpt(IOobject::NO_READ);
    }

    return io;
}


void Foam::radiation::radiationModel::initialise()
{
    if (radiation_)
    {
        solverFreq_ = max(1, getOrDefault<label>("solverFreq", 1));

        if (this->found("absorptionEmissionModel"))
        {
            absorptionEmission_.reset
            (
                absorptionEmissionModel::New(*this, mesh_).ptr()
            );
        }

        if (this->found("scatterModel"))
        {
            scatter_.reset(scatterModel::New(*this, mesh_).ptr());
        }

        if (this->found("sootModel"))
        {
            soot_.reset(sootModel::New(*this, mesh_).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::radiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(false),
    coeffs_(),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{}


Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(getOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    if (readOpt() == IOobject::NO_READ)
    {
        radiation_ = false;
    }

    initialise();
}


Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(getOrDefault("radiation", true)),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::radiationModel::read()
{
    if (regIOobject::read())
    {
        readEntry("radiation", radiation_);
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = getOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }

    return false;
}


void Foam::radiation::radiationModel::correct()
{
    if (!radiation_)
    {
        return;
    }

    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }

    if (soot_)
    {
        soot_->correct();
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::Sh
(
    const basicThermo& thermo,
    const volScalarField& he
) const
{
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ - 4.0*he/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::ST
(
    tmp<volScalarField> rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp.ref()
      - fvm::Sp(Rp()*pow3(T)/rhoCp.ref(), T)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::ST
(
    volScalarField& T
) const
{
    return
    (
        Ru()
      - fvm::Sp(Rp()*pow3(T), T)
    );
}


const Foam::radiation::absorptionEmissionModel&
Foam::radiation::radiationModel::absorptionEmission() const
{
    if (!absorptionEmission_)
    {
        FatalErrorInFunction
            << "Requested radiation absorptionEmission model, but model is "
            << "not activate" << abort(FatalError);
    }

    return *absorptionEmission_;
}


const Foam::radiation::sootModel&
Foam::radiation::radiationModel::soot() const
{
    if (!soot_)
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not activate" << abort(FatalError);
    }

    return *soot_;
}

/*
const Foam::radiation::transmissivityModel&
Foam::radiation::radiationModel::transmissivity() const
{
    if (!transmissivity_)
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not activate" << abort(FatalError);
    }

    return *transmissivity_;
}
*/

// ************************************************************************* //

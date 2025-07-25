/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "laminar.H"
#include "fvmSup.H"
#include "localEulerDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::laminar<ReactionThermo>::laminar
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    ChemistryCombustion<ReactionThermo>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    integrateReactionRate_
    (
        this->coeffs().getOrDefault("integrateReactionRate", true)
    )
{
    if (integrateReactionRate_)
    {
        Info<< "    using integrated reaction rate" << endl;
    }
    else
    {
        Info<< "    using instantaneous reaction rate" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::combustionModels::laminar<ReactionThermo>::~laminar()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminar<ReactionThermo>::tc() const
{
    return this->chemistryPtr_->tc();
}


template<class ReactionThermo>
void Foam::combustionModels::laminar<ReactionThermo>::correct()
{
    if (this->active())
    {
        if (integrateReactionRate_)
        {
            if (fv::localEulerDdt::enabled(this->mesh()))
            {
                const scalarField& rDeltaT =
                    fv::localEulerDdt::localRDeltaT(this->mesh());

                scalar maxTime;
                if (this->coeffs().readIfPresent("maxIntegrationTime", maxTime))
                {
                    this->chemistryPtr_->solve
                    (
                        min(1.0/rDeltaT, maxTime)()
                    );
                }
                else
                {
                    this->chemistryPtr_->solve((1.0/rDeltaT)());
                }
            }
            else
            {
                this->chemistryPtr_->solve(this->mesh().time().deltaTValue());
            }
        }
        else
        {
            this->chemistryPtr_->calculate();
        }
    }
}


template<class ReactionThermo>
Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::laminar<ReactionThermo>::R(volScalarField& Y) const
{
    auto tSu = tmp<fvScalarMatrix>::New(Y, dimMass/dimTime);
    auto& Su = tSu.ref();

    if (this->active())
    {
        const label specieI =
            this->thermo().composition().species().find(Y.member());

        Su += this->chemistryPtr_->RR(specieI);
    }

    return tSu;
}


template<class ReactionThermo>
Foam::tmp<Foam::volScalarField>
Foam::combustionModels::laminar<ReactionThermo>::Qdot() const
{
    auto tQdot = volScalarField::New
    (
        this->thermo().phaseScopedName(typeName, "Qdot"),
        IOobject::NO_REGISTER,
        this->mesh(),
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );

    if (this->active())
    {
        tQdot.ref() = this->chemistryPtr_->Qdot();
    }

    return tQdot;
}


template<class ReactionThermo>
bool Foam::combustionModels::laminar<ReactionThermo>::read()
{
    if (ChemistryCombustion<ReactionThermo>::read())
    {
        integrateReactionRate_ =
            this->coeffs().getOrDefault("integrateReactionRate", true);
        return true;
    }

    return false;
}


// ************************************************************************* //

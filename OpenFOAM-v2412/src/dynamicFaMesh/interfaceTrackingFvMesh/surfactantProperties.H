/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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

Class
    surfactantProperties

Description

SourceFiles
    surfactantProperties.H

\*---------------------------------------------------------------------------*/

#ifndef surfactantProperties_H
#define surfactantProperties_H

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class surfactantProperties Declaration
\*---------------------------------------------------------------------------*/

class surfactantProperties
{
    // Private data

        //- Surfactant concentration in the bulk of fluid
        dimensionedScalar bulkConc_;

        //- Saturated surfactant concentration on the free-surface
        dimensionedScalar saturatedConc_;

        //- Adsorption coefficient of surfactant
        dimensionedScalar adsorptionCoeff_;

        //- Desorption coefficient of surfactant
        dimensionedScalar desorptionCoeff_;

        //- Diffusion coefficient of surfactant in the bulk of fluid
        dimensionedScalar bulkDiffusion_;

        //- Diffusion coefficient of surfactant at the free-surface
        dimensionedScalar diffusion_;

        //- Temperature of surfactant at the free-surface
        dimensionedScalar T_;

        //- Universal gas constant
        dimensionedScalar R_;

        //- Equilibrium surfactant concentration at the free-surface
        dimensionedScalar equilibriumConc_;

        //- Is the surfactant soluble?
        Switch soluble_;


public:

    // Constructors

    explicit surfactantProperties(const dictionary& dict)
    :
        bulkConc_("bulkConc", dict),
        saturatedConc_("saturatedConc", dict),
        adsorptionCoeff_("adsorptionCoeff", dict),
        desorptionCoeff_("desorptionCoeff", dict),
        bulkDiffusion_("bulkDiffusion", dict),
        diffusion_("diffusion", dict),
        T_("temperature", dict),
        R_("R", dimGasConstant*dimMass/dimMoles, 8.3144),
        equilibriumConc_
        (
            saturatedConc_
           /(
                1.0
              + desorptionCoeff_
               /bulkConc_
            )
        ),
        soluble_(dict.get<bool>("soluble"))
    {}


    // Member Functions

        //- Return surfactant concentration in the bulk of fluid
        const dimensionedScalar& bulkConc() const
        {
            return bulkConc_;
        }

        //- Return saturated surfactant concentration at the free-surface
        const dimensionedScalar& saturatedConc() const
        {
            return saturatedConc_;
        }

        //- Return surfactant adsorption coefficient
        const dimensionedScalar& adsorptionCoeff() const
        {
            return adsorptionCoeff_;
        }

        //- Return surfactant desorption coefficient
        const dimensionedScalar& desorptionCoeff() const
        {
            return desorptionCoeff_;
        }

        //- Return diffusion coefficient of the surfactant in the bulk of fluid
        const dimensionedScalar& bulkDiffusion() const
        {
            return bulkDiffusion_;
        }

        //- Return diffusion coefficient of the surfactant at the free-surface
        const dimensionedScalar& diffusion() const
        {
            return diffusion_;
        }

        //- Return surfactant temeprature
        const dimensionedScalar& T() const
        {
            return T_;
        }

        //- Return universal gas constant
        const dimensionedScalar& R() const
        {
            return R_;
        }

        //- Return equilibrium surfactant concentration at the free-surface
        const dimensionedScalar& equilibriumConc() const
        {
            return equilibriumConc_;
        }

        //- Is the surfactant soluble
        Switch soluble() const
        {
            return soluble_;
        }

        //- Surface tension change due to presence of surfactants
        tmp<areaScalarField> dSigma
        (
            const areaScalarField& surfactConc
        ) const
        {
            return tmp<areaScalarField>::New
            (
                IOobject
                (
                    "dSigma",
                    surfactConc.mesh().time().timeName(),
                    surfactConc.mesh().thisDb(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                (
                    R()*T()*saturatedConc()
                  * log(1.0 - surfactConc/saturatedConc())
                )
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

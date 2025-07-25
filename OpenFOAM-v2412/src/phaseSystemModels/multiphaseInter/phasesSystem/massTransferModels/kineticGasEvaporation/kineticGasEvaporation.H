/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::meltingEvaporationModels::kineticGasEvaporation

Description
    Considering the Hertz Knudsen formula, which gives the
    evaporation-condensation flux based on the kinetic theory for flat
    interface:

    \f[
        Flux = C \sqrt{\frac{M}{2 \pi R T_{activate}}}(p - p_{sat})
    \f]

        where:
        \vartable
            Flux    | mass flux rate [kg/s/m2]
            M       | molecular weight
            T_{activate} | saturation temperature
            C       | accommodation coefficient
            R       | universal gas constant
            p_{sat} | saturation pressure
            p       | vapor partial pressure
        \endvartable

    The Clapeyron-Clausius equation relates the pressure to the temperature
    for the saturation condition:

    \f[
        \frac{dp}{dT} = - \frac{L}{T (\nu_v - \nu_l)}
    \f]

        where:
        \vartable
            L       | latent heat
            \nu_v   | inverse of the vapor density
            \nu_l   | inverse of the liquid density
         \endvartable


    Using the above relations:

    \f[
        Flux =
           2 \frac{C}{2 - C}
           \sqrt{\frac{M}{2 \pi R {T_activate}^3}} L (\rho_{v})
           (T - T_{activate})
    \f]

    This assumes liquid and vapour are in equilibrium, then the accommodation
    coefficient are equivalent for the interface. This relation is known as the
    Hertz-Knudsen-Schrage.

    Based on the reference:
    - Van P. Carey, Liquid-Vapor Phase Change Phenomena, ISBN 0-89116836,
      1992, pp. 112-121.


Usage

    Example usage:
    \verbatim
    massTransferModel
    (
        (liquid to gas)
        {
            type                kineticGasEvaporation;
            species             vapour.gas;
            C                   0.1;
            isoAlpha            0.1;
            Tactivate           373;
        }
    );
    \endverbatim

    where:
    \table
        Property    | Description             | Required    | Default value
        C           | Coefficient (C > 0 for evaporation, C < 0 for
            condensation) | yes
        includeVolChange    | Volumen change  | no          | yes
        isoAlpha    | iso-alpha for interface | no          | 0.5
        Tactivate   | Saturation temperature  | yes
        species     | Specie name on the other phase | no   | none
    \endtable

SourceFiles
    kineticGasEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_kineticGasEvaporation_H
#define meltingEvaporationModels_kineticGasEvaporation_H

#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

// Forward Declarations
class phasePair;

namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                    Class kineticGasEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class kineticGasEvaporation
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private Data

        //- Evaporation coefficient
        dimensionedScalar C_;

        //- Activation temperature
        const dimensionedScalar Tactivate_;

        //- Molar weight of the vapour in the continuous phase
        dimensionedScalar Mv_;

        //- Interface area
        volScalarField interfaceArea_;

        //- Heat transfer coefficient
        volScalarField htc_;

        //- Mass source
        volScalarField mDotc_;

        //- Interface Iso-value
        scalar isoAlpha_;


    // Private Member Functions

        //- Update interface
        void updateInterface(const volScalarField& T);


public:

    //- Runtime type information
    TypeName("kineticGasEvaporation");


    // Constructors

        //- Construct from components
        kineticGasEvaporation
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~kineticGasEvaporation() = default;


    // Member Functions

        //- Explicit total mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            const volScalarField& field
        );

        //- Implicit mass transfer coefficient
        virtual tmp<volScalarField> KSp
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> KSu
        (
            label modelVariable,
            const volScalarField& field
        );

        //- Return Tactivate
        virtual const dimensionedScalar& Tactivate() const noexcept
        {
            return Tactivate_;
        }

        //- Add/subtract alpha*div(U) as a source term
        //- for alpha, substituting div(U) = mDot(1/rho1 - 1/rho2)
        virtual bool includeDivU() const noexcept
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kineticGasEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //

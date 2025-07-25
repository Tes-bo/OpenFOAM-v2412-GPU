/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::meltingEvaporationModels::diffusionGasEvaporation

Description
    Gas diffusion based evaporation/condensation mass transfer model.

    THE evaporation rate is given by:

    \f[
        \hat{m_i} =
            - C \rho_{g} D_{v,i}
            \frac
            {
                \frac{dY_{v,i}}{dn}
            }
            {
                1 - \; \sum_{1}^{Nv}Y_{g,j}
            }
    \f]

    where:
    \vartable
        \hat{m_i}           | mass flux rate [kg/s/m2]
        \rho_{g}            | gas phase density
        D_{v,i}             | diffusion coefficient
        C                   | model coefficient
        \frac{dY_{v,i}}{dn} | normal derivative of evaporated component
        Y_{g,j}             | mass fraction at the surface
    \endvartable

Usage
    Example usage:
    \verbatim
    massTransferModel
    (
        (liquid to gas)
        {
            // Mandatory entries
            type                diffusionGasEvaporation;
            species             vapour.gas;
            C                   1;
            saturationPressure
            {
                type    Antoine;
                A       3.55;
                B       643;
                C       -198;
            }

            // Optional entries
            isoAlpha            <scalar>;
            Tactivate           <scalar>;

            // Inherited entries
            ...
        }
    );
    \endverbatim

    where:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      type         | Type name: diffusionGasEvaporation | word | yes  | -
      saturationPressure | Saturation model             | dict | yes  | -
      isoAlpha    | Iso-alpha for interface            | scalar | no  | 0.5
      C           | Model coefficient                  | scalar | yes | -
      Tactivate   | Saturation temperature             | scalar | no  | 0
    \endtable

    The inherited entries are elaborated in:
      - \link InterfaceCompositionModel.H \endlink
      - \link saturationModel.H \endlink

SourceFiles
    diffusionGasEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_meltingEvaporationModels_diffusionGasEvaporation_H
#define Foam_meltingEvaporationModels_diffusionGasEvaporation_H

#include "InterfaceCompositionModel.H"
#include "saturationModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

// Forward Declarations
class phasePair;

namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                    Class diffusionGasEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class diffusionGasEvaporation
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private Data

        //- Saturation pressure model
        autoPtr<saturationModel> saturationModelPtr_;

        //- Interface Iso-value
        scalar isoAlpha_;

        //- Evaporation coefficient
        dimensionedScalar C_;

        //- Activation temperature
        const dimensionedScalar Tactivate_;

        //- Interface area
        volScalarField interfaceArea_;

        //- Mass source
        volScalarField mDotc_;


    // Private Member Functions

        //- Update interface
        void updateInterface(const volScalarField& T);


public:

    //- Runtime type information
    TypeName("diffusionGasEvaporation");


    // Constructors

        //- Construct from components
        diffusionGasEvaporation
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~diffusionGasEvaporation() = default;


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
#   include "diffusionGasEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //

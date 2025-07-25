/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::meltingEvaporationModels::interfaceOxideRate

Description
    The \c interfaceOxideRate is a simple model to calculate the formation
    rate of oxide inclusions (\c mDotOxide), which is a local function of
    the volume fraction of reducing agent (\c alpha), temperature (\c T)
    and oxide-inclusion density (\c chi).

    The oxide-inclusion formation rate is modelled as follows:

    \f[
        S_\chi = C_\chi S_\alpha S_T S_\rho
    \f]

    with

    \f[
        S_\alpha = \alpha (1 - \alpha)
    \f]

    \f[
        S_T =
        \exp{
            1 - \frac{1}{max(T - T_{solidus}/(T_liquidus - T_solidus), 1e-6)}
        }
    \f]

    \f[
        S_\rho =
            max \left(\frac{\chi_{crit} - \chi_{curr}}{\chi_{curr}}, 0\right)
    \f]

    where
    \vartable
      S_\chi   | Oxide-inclusion formation rate          [kg/(m^3 s)]
      C_\chi   | Oxide-inclusion formation rate constant [kg/(m^3 s)]
      S_\alpha | Formation factor due to volume fraction of reducing agent [-]
      S_T      | Formation factor due to temperature     [-]
      S_\rho   | Formation factor due to oxide-inclusion density [-]
      \alpha   | Volume fraction of reducing agent       [-]
      T        | Local temperature                       [K]
      T_{solidus}  | Solidus temperature of reducing agent  [K]
      T_{liquidus} | Liquidus temperature of reducing agent [K]
      \chi_{crit}  | Critical oxide-inclusion density    [kg/m^3]
      \chi_{curr}  | Current oxide-inclusion density     [kg/m^3]
    \endvartable

    References:
    \verbatim
        Oxide-inclusion model (tag:CSC):
            Cao, L., Sun, F., Chen, T., Tang, Y., & Liao, D. (2018).
            Quantitative prediction of oxide inclusion defects inside
            the casting and on the walls during cast-filling processes.
            International Journal of Heat and Mass Transfer, 119, 614-623.
            DOI:10.1016/j.ijheatmasstransfer.2017.11.127
    \endverbatim

Usage
    Minimal example by using \c constant/phaseProperties.massTransferModel:
    \verbatim
    massTransferModel
    (
        (liquid to oxide)
        {
            type                interfaceOxideRate;
            C                   <scalar>;
            Tliquidus           <scalar>;
            Tsolidus            <scalar>;
            oxideCrit           <scalar>;
            isoAlpha            <scalar>;
        }
    );
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                             | Type | Reqd | Deflt
      type      | Type name: interfaceOxideRate           | word   | yes | -
      C         | Oxide-inclusion formation rate constant | scalar | yes | -
      Tliquidus | Liquidus temperature of reducing agent  | scalar | yes | -
      Tsolidus  | Solidus temperature of reducing agent   | scalar | yes | -
      oxideCrit | Critical oxide-inclusion density        | scalar | yes | -
      isoAlpha  | Location of the source                  | scalar | no  | 0.5
    \endtable

Note
  - \c oxideCrit should be determined experimentally (CSC:p. 616).
  - \c C should be determined by practical production (CSC:p. 616).

SourceFiles
    interfaceOxideRate.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_interfaceOxideRate_H
#define meltingEvaporationModels_interfaceOxideRate_H

#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

// Forward Declarations
class phasePair;

namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                    Class interfaceOxideRate Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class interfaceOxideRate
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private Data

        //- Oxide-inclusion formation rate constant
        const dimensionedScalar C_;

        //- Liquidus temperature of reducing agent (e.g. a casting metal)
        const dimensionedScalar Tliquidus_;

        //- Solidus temperature of reducing agent (e.g. a casting metal)
        const dimensionedScalar Tsolidus_;

        //- Critical oxide-inclusion density
        const dimensionedScalar oxideCrit_;

        //- Oxide-inclusion formation rate
        volScalarField mDotOxide_;

        //- Interface Iso-value
        scalar isoAlpha_;


public:

    //- Runtime type information
    TypeName("interfaceOxideRate");


    // Constructors

        //- Construct from components
        interfaceOxideRate
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceOxideRate() = default;


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

        //- Return solidus temperature of reducing agent
        virtual const dimensionedScalar& Tactivate() const noexcept
        {
            return Tsolidus_;
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
#   include "interfaceOxideRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //

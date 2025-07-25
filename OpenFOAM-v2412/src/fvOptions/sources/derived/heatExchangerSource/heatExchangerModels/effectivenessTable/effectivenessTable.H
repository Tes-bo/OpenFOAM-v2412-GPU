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
    Foam::heatExchangerModels::effectivenessTable

Description
    A heat exchanger model where heat exchange
    is calculated via an effectiveness table.

    The total heat exchange source is given by:
    \f[
        Q_t = e(\phi, \dot{m}_2) (T_2 - T_1) \phi C_p
    \f]

    where:
    \vartable
        Q_t       | Total heat exchange source [J/s]
        e(\phi,\dot{m}_2) | Effectivenes table [-]
        \phi      | Net mass flux entering heat exchanger [kg/s]
        \dot{m}_2 | Secondary flow mass flow rate [kg/s]
        T_1       | Primary flow inlet temperature [K]
        T_2       | Secondary flow inlet temperature [K]
        C_p       | Primary flow specific heat capacity [J/kg/K]
    \endvartable

    The distribution inside the heat exchanger is given by:
    \f[
        Q_c = \frac{V_c |U_c| (T_c - T_{ref})}{\sum(V_c |U_c| (T_c - T_{ref}))}
    \f]

    where:
    \vartable
        Q_c     | Source for cell
        V_c     | Volume of the cell [m^3]
        U_c     | Local cell velocity [m/s]
        T_c     | Local cell temperature [K]
        T_{ref} | Min or max(T) in cell zone depending on the sign of Qt [K]
    \endvartable

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <name>
    {
        // Inherited entries
        ...

        // Mandatory entries
        model                    effectivenessTable;
        secondaryMassFlowRate    <scalar>;
        secondaryInletT          <scalar>;
        file                     "<word>";
        outOfBounds              <word>;

        // Conditional optional entries

            // when the total heat exchange is calculated with primary inlet T
            primaryInletT           <scalar>;

            // when the total heat exchange is calculated with a given target
            targetQdot              <scalar>;
            targetQdotCalcInterval  <label>;
            targetQdotRelax         <scalar>;

            // when secondary outlet temperature is requested
            secondaryCp             <Function1<scalar>>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type  | Reqd | Deflt
      model     | Model name:effectivenessTable     | word  | yes  | -
      secondaryMassFlowRate | Secondary flow mass rate [kg/s]     <!--
                -->                                 | scalar | yes | -
      secondaryInletT | Secondary flow inlet temperature [K]      <!--
                -->                                 | scalar | yes | -
      file      | 2D effectiveness table = function of primary    <!--
                --> and secondary mass flow rates [kg/s] | file | yes | -
      primaryInletT | Primary flow temperature at the heat exchanger inlet <!--
                -->                                 | scalar | no  | -
      targetQdot |  Target heat rejection           | scalar | no  | -
      targetQdotCalcInterval | Target heat rejection calculation interval <!--
                -->                                 | label  | no  | -
      targetQdotRelax | Target heat rejection temperature <!--
                --> under-relaxation coefficient    | scalar | no  | -
      secondaryCp | Secondary flow specific heat capacity <!--
                -->                    | Function1\<scalar\> | no  | -
    \endtable

    The inherited entries are elaborated in:
      - \link heatExchangerSource.H \endlink
      - \link interpolation2DTable.H \endlink

    The effectiveness table is described in terms of the primary and secondary
    mass flow rates.  For example, the table:

    \verbatim
                           secondary MFR
                       |  0.1   0.2   0.3
                  -----+-----------------
                  0.02 |   A     B     C
     primary MFR  0.04 |   D     E     F
                  0.06 |   G     H     I
    \endverbatim

    is specified by the following:

    \verbatim
        (
            (
                0.02
                (
                    (0.1    A)
                    (0.2    B)
                    (0.3    C)
                )
            )
            (
                0.04
                (
                    (0.1    D)
                    (0.2    E)
                    (0.3    F)
                )
            )
            (
                0.06
                (
                    (0.1    G)
                    (0.2    H)
                    (0.3    I)
                )
            )
        );
    \endverbatim

Note
  - Primary flow indicates the CFD flow region and
    secondary flow the non-CFD-model region.
  - The table with name \c file should have the same units as the
    secondary mass flow rate and kg/s for \c phi.
  - \c faceZone is the faces at the inlet of the \c cellZone, it needs to be
    created with flip map flags. It is used to integrate the net mass flow
    rate into the heat exchanger.
  - \c primaryInletT sets the primary inlet temperature.  If not set, the
    flux-averaged temperature is used.

SourceFiles
    effectivenessTable.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_heatExchangerModels_effectivenessTable_H
#define Foam_heatExchangerModels_effectivenessTable_H

#include "heatExchangerModel.H"
#include "interpolation2DTable.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatExchangerModels
{

/*---------------------------------------------------------------------------*\
                              Class effectivenessTable Declaration
\*---------------------------------------------------------------------------*/

class effectivenessTable
:
    public heatExchangerModel
{
    // Private Data

        //- Flag to use a user-specified primary flow inlet temperature
        bool userPrimaryInletT_;

        //- Flag to use target heat rejection
        bool targetQdotActive_;

        //- Secondary flow specific heat capacity [J/kg/K]
        autoPtr<Function1<scalar>> secondaryCpPtr_;

        //- 2D effectiveness table = function of primary and secondary
        //- mass flow rates [kg/s]
        autoPtr<interpolation2DTable<scalar>> eTable_;

        //- Target heat rejection calculation interval
        label targetQdotCalcInterval_;

        //- Secondary flow mass rate [kg/s]
        scalar secondaryMassFlowRate_;

        //- Secondary flow inlet temperature [K]
        scalar secondaryInletT_;

        //- Primary flow temperature at the heat exchanger inlet [K]
        scalar primaryInletT_;

        //- Target heat rejection
        scalar targetQdot_;

        //- Target heat rejection temperature under-relaxation coefficient
        scalar targetQdotRelax_;

        //- Net mass flux [kg/s]
        scalar sumPhi_;

        //- Total heat exchange [W]
        scalar Qt_;

        //- Reference temperature [K]
        scalar Tref_;

        //- Effectiveness of the heat exchanger [-]
        scalar effectiveness_;


    // Private Member Functions

    // I-O

        //- Write file header information
        void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("effectivenessTable");


    // Constructors

        //- Construct from components
        effectivenessTable
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& coeffs
        );

        //- No copy construct
        effectivenessTable(const effectivenessTable&) = delete;

        //- No copy assignment
        void operator=(const effectivenessTable&) = delete;


    //- Destructor
    virtual ~effectivenessTable() = default;


    // Member Functions

    // Evaluation

        //- Initialise data members of the model
        virtual void initialise();

        //- Return energy density per unit length [J/m3/m]
        virtual tmp<scalarField> energyDensity(const labelList& cells);


    // I-O

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict);

        //- Write data to stream and files
        virtual void write(const bool log);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatExchangerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

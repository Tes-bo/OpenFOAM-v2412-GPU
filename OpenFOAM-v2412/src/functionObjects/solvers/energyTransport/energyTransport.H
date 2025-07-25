/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::functionObjects::energyTransport

Group
    grpSolversFunctionObjects

Description
    Computes the simplified energy transport equation in single-phase or
    two-phase flow, considering incompressible cases:

    \f[
        \frac{\partial \rho \, C_p \, T}{\partial t}
      + \nabla \cdot \left(\rho \, C_p \, \phi \, T \right)
      - \nabla \cdot \left(\rho \, C_p \, \phi \right) \, T
      - \nabla \cdot \left(\kappa_{eff} \, \nabla T \right)
      = S_T
    \f]

    where:
    \vartable
      T       | Scalar field
      \rho    | (Generic) Fluid density which is unity when not specified
      C_p     | Specific heat capacity at constant pressure
      \phi    | (Generic) Flux field
      \kappa_{eff}  | Effective thermal conductivity
      S_T     | Scalar field source term
    \endvartable

Usage
    Minimal example in \c system/controlDict.functions:
    \verbatim
    energyTransport1
    {
        // Mandatory entries
        type            energyTransport;
        libs            (solverFunctionObjects);

        // Optional entries
        field           <word>;
        phi             <word>;
        rho             <word>;
        Cp              <scalar>;
        kappa           <scalar>;
        rhoInf          <scalar>;
        Prt             <scalar>;
        schemesField    <word>;
        tolerance       <scalar>;
        nCorr           <int>;
        fvOptions       <dict>;
        phaseThermos    <dict>;

        // Inherited entries
        ...
    }
    \endverbatim

    where:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: energyTransport          | word | yes  | -
      libs      | Library name: solverFunctionObjects | word | yes  | -
      field     | Name of the passive-scalar field    | word | no   | s
      phi       | Name of flux field                  | word | no   | phi
      rho       | Name of density field               | word | no   | rho
      Cp        | Specific heat capacity at constant pressure | scalar | no | 0
      kappa     | Thermal conductivity                | scalar | no | 0
      rhoInf    | Fluid density                       | scalar | no | 0
      Prt       | Turbulent Prandtl number            | scalar | no | 1
      schemesField | Name of field to specify schemes | word   | no | field
      tolerance | Outer-loop initial-residual tolerance | scalar | no | 1
      nCorr     | Number of outer-loop correctors     | int    | no | 0
      fvOptions | List of finite-volume options       | dict   | no | -
      phaseThermos | Dictionary for multi-phase thermo | dict  | no | null
    \endtable

    The inherited entries are elaborated in:
      - \link fvMeshFunctionObject.H \endlink
      - \link fvOption.H \endlink

    An example of function object specification to solve a energy transport
    equation for a single phase flow plus a source term:
    \verbatim
    energyTransport1
    {
        // Mandatory entries
        type            energyTransport;
        libs            (solverFunctionObjects);

        // Optional entries
        field           T;
        phi             phi;
        Cp              Cp    [J/kg/K]  1e3;
        kappa           kappa [W/m/K]   0.0257;
        rhoInf          rho   [kg/m^3]  1.2;
        fvOptions
        {
            viscousDissipation
            {
                type            viscousDissipation;
                enabled         true;

                viscousDissipationCoeffs
                {
                    fields          (T);
                    rhoInf          $....rhoInf;
                }
            }
        }

        // Inherited entries
        enabled         true;
        writeControl    writeTime;
        writeInterval   1;
    }
    \endverbatim

    An example of function object specification to solve a energy transport
    equation for a multiphase phase flow plus a source term:
    \verbatim
    energyTransport1
    {
        // Mandatory entries
        type            energyTransport;
        libs            (solverFunctionObjects);

        // Optional entries
        field           T;
        rho             rho;
        phi             rhoPhi;

        // Thermal properties of the phases
        phaseThermos
        {
            alpha.air
            {
                Cp          1e3;
                kappa       0.0243;
            }
            alpha.mercury
            {
                Cp          140;
                kappa       8.2;
            }
            alpha.oil
            {
                Cp          2e3;
                kappa       0.2;
            }
            alpha.water
            {
                Cp          4e3;
                kappa       0.6;
            }
        }

        fvOptions
        {
            viscousDissipation
            {
                type            viscousDissipation;
                enabled         true;

                viscousDissipationCoeffs
                {
                    fields          (T);
                    rho             rho;
                }
            }
        }

        // Inherited entries
        enabled         true;
        writeControl    writeTime;
        writeInterval   1;
    }
    \endverbatim

Note
  - The field name must be temperature and its boundary conditions
    specified in the time directory.
  - The turbulence model should be incompressible.

SourceFiles
    energyTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_energyTransport_H
#define functionObjects_energyTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class energyTransport Declaration
\*---------------------------------------------------------------------------*/

class energyTransport
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Volumetric heat capacity field [J/m^3/K]
        volScalarField rhoCp_;

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- Dictionary for multiphase thermos
        dictionary multiphaseThermo_;

        //- List of phase names
        wordList phaseNames_;

        //- List of phase specific heat capacities at constant pressure
        PtrList<dimensionedScalar> Cps_;

        //- List of phase thermal diffusivity for temperature [J/m/s/K]
        PtrList<dimensionedScalar> kappas_;

        //- Unallocated phase list
        UPtrList<volScalarField> phases_;

        //- Specific heat capacity at constant pressure for single phase flows
        dimensionedScalar Cp_;

        //- Thermal diffusivity for temperature for single phase flows
        dimensionedScalar kappa_;

        //- Density for single phase flows
        dimensionedScalar rho_;

        //- Turbulent Prandt number
        dimensionedScalar Prt_;

        //- Name of the transport field
        word fieldName_;

        //- Name of field whose schemes are used
        word schemesField_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Outer-loop initial-residual tolerance
        scalar tol_;

        //- Number of corrector iterations
        int nCorr_;


    // Private Member Functions

        //- Return reference to registered transported field
        volScalarField& transportedField();

        //- Return the diffusivity field
        tmp<volScalarField> kappaEff() const;

        //- Return the density field, rho
        tmp<volScalarField> rho() const;

        //- Return the specific heat capacity at constant pressure field, Cp
        tmp<volScalarField> Cp() const;

        //- Return the thermal diffusivity field
        tmp<volScalarField> kappa() const;


public:

    //- Runtime type information
    TypeName("energyTransport");


    // Constructors

        //- Construct from Time and dictionary
        energyTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~energyTransport() = default;


    // Member Functions

        //- Read the energyTransport data
        virtual bool read(const dictionary&);

        //- Calculate the energyTransport
        virtual bool execute();

        //- Do nothing.
        //  The volScalarField is registered and written automatically
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

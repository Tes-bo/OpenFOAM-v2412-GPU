/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::functionObjects::electricPotential

Group
    grpSolversFunctionObjects

Description
    Computes the steady-state equation of charge conservation to obtain
    the electric potential by strictly assuming a quasi-static electrostatic
    field for single-phase and multiphase applications.

    The steady-state equation of the charge conservation:

    \f[
        \nabla \cdot \left( \sigma \nabla V \right) = 0
    \f]

    where
    \vartable
      V       | Electric potential                [volt = kg m^2/(A s^3)]
      \sigma  | Isotropic conductivity of mixture [S/m = A^2 s^3/(kg m^3)]
    \endvartable

    Optionally, electric field, current density and free-charge
    density fields can be written out by using the following equations:

    \f[
        \vec{E} = - \nabla V
    \f]

    \f[
        \vec{J} = \sigma \vec{E} = - \sigma \nabla V
    \f]

    \f[
        \rho_E = \nabla \cdot \left(\epsilon_m \vec{E} \right)
               = \nabla \cdot \left(\epsilon_0 \epsilon_r \vec{E} \right)
    \f]

    where
    \vartable
      \vec{E}    | Electric field                       [m kg/(s^3 A)]
      \vec{J}    | Current density                      [A/m^2]
      \rho_E     | Volume charge density                [C/m^3 = A s/m^3]
      \epsilon_m | Isotropic permittivity of mixture    [F/m = A^2 s^4/(kg m^3)]
      \epsilon_0 | Isotropic vacuum permittivity        [F/m = A^2 s^4/(kg m^3)]
      \epsilon_r | Isotropic relative permittivity of mixture [-]
    \endvartable

    For multiphase applications, \c sigma and \c epsilonr are blended
    (to consider their interface values) by using the simple weighted
    arithmetic mean interpolation, for example:

    \f[
        \sigma = \alpha_1 \sigma_1 + \alpha_2 \sigma_2
               = \alpha_1 \sigma_1 + (1 - \alpha_1) \sigma_2
    \f]

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    electricPotential1
    {
        // Mandatory entries
        type            electricPotential;
        libs            (solverFunctionObjects);

        // Conditional entries

            // Option-1: single-phase
            sigma       <scalar>;
            epsilonr    <scalar>;

            // Option-2: multiphase
            phases
            {
                alpha.air
                {
                    sigma       <scalar>;
                    epsilonr    <scalar>;
                }
                alpha.water
                {
                    sigma       <scalar>;
                    epsilonr    <scalar>;
                }
                alpha.mercury
                {
                    sigma       <scalar>;
                    epsilonr    <scalar>;
                }
                ...
            }

        // Optional entries
        nCorr                 <int>;
        writeDerivedFields    <bool>;
        V                     <word>;
        electricField         <bool>;
        E                     <word>;
        fvOptions             <dict>;
        tolerance             <scalar>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: electricPotential        | word | yes  | -
      libs      | Library name: solverFunctionObjects | word | yes  | -
      sigma     | Isotropic electrical conductivity of phase | scalar | yes | -
      epsilonr  | Isotropic relative permittivity of phase   | scalar | no  | -
      nCorr     | Number of corrector iterations      | int  | no   | 1
      writeDerivedFields | Flag to write extra fields | bool | no   | false
      V    | Name of electric potential field | word | no | electricPotential:V
      electricField | Flag to calculate electric field  | bool | no | false
      E    | Name of electric field           | word | no | electricPotential:E
      fvOptions | List of finite-volume options       | dict | no   | -
      tolerance | Outer-loop initial-residual tolerance | scalar | no | 1
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link fvOption.H \endlink

    Fields written out when the \c writeDerivedFields entry is \c true:
    \table
      Operand         | Type             | Location
      Current density | volVectorField   | \<time\>/electricPotential:J
      Charge density  | volScalarField   | \<time\>/electricPotential:rho
    \endtable

Note
  - Only constraint-type finite-volume options can be used.

SourceFiles
    electricPotential.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_electricPotential_H
#define functionObjects_electricPotential_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class electricPotential Declaration
\*---------------------------------------------------------------------------*/

class electricPotential
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Dictionary of phase data
        dictionary phasesDict_;

        //- List of phase names
        wordList phaseNames_;

        //- Unallocated list of phase fields
        UPtrList<volScalarField> phases_;

        //- List of isotropic electrical conductivity of phases
        PtrList<dimensionedScalar> sigmas_;

        //- Isotropic electrical conductivity of a single phase
        dimensionedScalar sigma_;

        //- List of isotropic relative permittivity of phases
        PtrList<dimensionedScalar> epsilonrs_;

        //- Isotropic relative permittivity of a single phase
        dimensionedScalar epsilonr_;

        //- Name of electric potential field
        word Vname_;

        //- Name of electric field
        word Ename_;

        //- Run-time selectable finite volume options
        fv::optionList fvOptions_;

        //- Outer-loop initial-residual tolerance
        scalar tol_;

        //- Number of corrector iterations
        int nCorr_;

        //- Flag to write derived fields of
        //- electric field, current density and free-charge density
        bool writeDerivedFields_;

        //- Flag to calculate electric field
        bool electricField_;


    // Private Member Functions

        //- Return requested field from the object registry
        //- or read+register the field to the object registry
        volScalarField& getOrReadField(const word& fieldName) const;


        //- Return the isotropic electrical conductivity field of the mixture
        tmp<volScalarField> sigma() const;

        //- Return the isotropic permittivity field of the mixture
        tmp<volScalarField> epsilonm() const;


        //- No copy construct
        electricPotential(const electricPotential&) = delete;

        //- No copy assignment
        void operator=(const electricPotential&) = delete;


public:

    //- Runtime type information
    TypeName("electricPotential");


    // Constructors

        //- Construct from Time and dictionary
        electricPotential
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~electricPotential() = default;


    // Member Functions

        //- Read the function object data
        virtual bool read(const dictionary& dict);

        //- Calculate the function object
        virtual bool execute();

        //- Write the function object output
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::functionObjects::scalarTransport

Group
    grpSolversFunctionObjects

Description
    Computes the transport equation for a passive scalar in single-phase or
    two-phase flow, considering both incompressible and compressible cases:

    \f[
        \frac{\partial \rho \, T}{\partial t}
      + \nabla \cdot \left( \phi_\alpha \, T \right)
      - \nabla \cdot (D_T \, \nabla T)
      = \alpha \, S_T
    \f]

    where:
    \vartable
      T       | Passive scalar field
      \rho    | (Generic) Fluid density which is unity when not specified
      \phi_\alpha | (Generic) Flux field
      \alpha  | Phase fraction which is unity for single-phase flows
      D_T     | Diffusivity representing the diffusive transport of T
      S_T     | Passive-scalar field source term
    \endvartable

Usage
    Minimal example in \c system/controlDict.functions:
    \verbatim
    scalarTransport1
    {
        // Mandatory entries
        type            scalarTransport;
        libs            (solverFunctionObjects);

        // Optional entries
        field           <word>;
        phi             <word>;
        rho             <word>;
        nut             <word>;
        phase           <word>;
        phasePhiCompressed  <word>;
        schemesField    <word>;
        bounded01       <bool>;
        D               <scalar>;
        alphaD          <scalar>;
        alphaDt         <scalar>;
        tolerance       <scalar>;
        nCorr           <int>;
        resetOnStartUp  <bool>;
        fvOptions       <dict>;

        // Inherited entries
        ...
    }

    where:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: scalarTransport          | word | yes  | -
      libs      | Library name: solverFunctionObjects | word | yes  | -
      field     | Name of the passive-scalar field    | word | no   | s
      phi       | Name of flux field                  | word | no   | phi
      rho       | Name of density field               | word | no   | rho
      nut       | Name of the turbulence viscosity    | word | no   | none
      phase     | Name of the phase                   | word | no   | none
      phasePhiCompressed | Name of compressed VOF flux | word | no | alphaPhiUn
      schemesField | Name of field to specify schemes | word  | no | field
      bounded01 | Bounds scalar between 0-1 for multiphase  | bool | no | true
      D         | Diffusion coefficient               | scalar | no | -
      alphaD    | Laminar diffusivity coefficient     | scalar | no | 1
      alphaDt   | Turbulent diffusivity coefficient   | scalar | no | 1
      tolerance | Outer-loop initial-residual tolerance | scalar | no | 1
      nCorr     | Number of outer-loop correctors     | int    | no | 0
      resetOnStartUp | Flag to reset field to zero on start-up | bool | no | no
      fvOptions | List of finite-volume options       | dict   | no | -
    \endtable

    The inherited entries are elaborated in:
      - \link fvMeshFunctionObject.H \endlink
      - \link fvOption.H \endlink

    An example of function object specification to solve a residence time
    in a two-phase flow:
    \verbatim
    scalarTransport1
    {
        // Mandatory entries
        type            scalarTransport;
        libs            (solverFunctionObjects);

        // Optional entries
        field           s;
        bounded01       false;
        phase           alpha.water;
        tolerance       1e-5;
        resetOnStartUp  false;
        fvOptions
        {
            unitySource
            {
                type            scalarSemiImplicitSource;
                enabled         true;

                selectionMode   all;
                volumeMode      specific;

                sources
                {
                    s           (1 0);
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
  - To use the same numerical schemes as another field,
  set the \c schemesField entry.
  - The diffusivity can be set manually using the \c D entry, obtained
  from the turbulence model or specified as `nut`.
  - Alternatively, if a turbulence model is available, turbulent diffusivity
  can be constructed from the laminar and turbulent viscosities using the
  optional diffusivity coefficients \c alphaD and \c alphaDt
  (which default to 1):

    \f[
        D = \alpha_D \, \nu + \alpha_{Dt} \, \nu_t
    \f]

SourceFiles
    scalarTransport.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_scalarTransport_H
#define functionObjects_scalarTransport_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "fvOptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class scalarTransport Declaration
\*---------------------------------------------------------------------------*/

class scalarTransport
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Run-time selectable finite volume options, e.g. sources, constraints
        fv::optionList fvOptions_;

        //- Name of the transport field.
        word fieldName_;

        //- Name of field whose schemes are used
        word schemesField_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Name of turbulent viscosity field
        word nutName_;

        //- Name of phase field
        word phaseName_;

        //- Name of phase field compressed flux
        word phasePhiCompressedName_;

        //- Diffusion coefficient
        scalar D_;

        //- Laminar diffusion coefficient
        scalar alphaD_;

        //- Turbulent diffusion coefficient
        scalar alphaDt_;

        //- Outer-loop initial-residual tolerance
        scalar tol_;

        //- Number of corrector iterations
        int nCorr_;

        //- Flag to reset the scalar to zero on start-up
        bool resetOnStartUp_;

        //- Flag to indicate whether a constant, uniform D_ is specified
        bool constantD_;

        //- Bound scalar between 0-1 using MULES for multiphase case
        bool bounded01_;


    // Private Member Functions

        //- Return reference to registered transported field
        volScalarField& transportedField();

        //- Return the diffusivity field
        tmp<volScalarField> D
        (
            const volScalarField& s,
            const surfaceScalarField& phi
        ) const;


public:

    //- Runtime type information
    TypeName("scalarTransport");


    // Constructors

        //- Construct from Time and dictionary
        scalarTransport
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~scalarTransport() = default;


    // Member Functions

        //- Read the scalarTransport data
        virtual bool read(const dictionary&);

        //- Calculate the scalarTransport
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

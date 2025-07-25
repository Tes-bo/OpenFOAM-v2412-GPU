/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd
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
    Foam::fv::viscousDissipation

Group
    grpFvOptionsSources

Description
    Applies a viscous dissipation energy source to the energy equation.

    Sources applied to either of the below, if exist:
    \verbatim
      e         | Internal energy                            [m2/s2]
      h         | Enthalphy                                  [m2/s2]
      field     | User-specified field
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity                                   [m/s]
      rho       | Fluid density (for compressible flows)     [kg/m3]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    viscousDissipation1
    {
        // Mandatory entries (unmodifiable)
        type                 viscousDissipation;

        // Optional entries (unmodifiable)
        fields           (<fieldName>);
        U                <UName>;
        rho              <rhoName>;

        // Conditional optional entries (unmodifiable)

            // when rho=rhoInf
            rhoInf       1000;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: viscousDissipation     | word   | yes   | -
      fields    | Name of operand (energy) field    | word   | no    | -
      U         | Name of operand velocity field    | word   | no    | U
      rho       | Name of operand density field     | word   | no    | none
      rhoInf    | Value of fluid density            | scalar | cndtnl | 0
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

Note
  - The entry \c fields accepts only a single field name.
  - The sources are applied onto the entire computational domain.
  - \c rho=rhoInf sets a zero-density for single phase flows.
  - Not applicable to multiphase flows.

SourceFiles
    viscousDissipation.C

\*---------------------------------------------------------------------------*/

#ifndef viscousDissipation_H
#define viscousDissipation_H

#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class viscousDissipation Declaration
\*---------------------------------------------------------------------------*/

class viscousDissipation
:
    public fv::option
{
    // Private Data

        //- Name of velocity field
        word UName_;

        //- Name of the rho field for incompressible solvers
        word rhoName_;

        //- Density for single phase flows
        dimensionedScalar rho_;


     // Private Member Functions

        //- Return the viscosity field
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return rho field
        tmp<volScalarField> rho() const;


public:

    //- Runtime type information
    TypeName("viscousDissipation");


    // Constructors

        //- Construct from explicit source name and mesh
        viscousDissipation
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        viscousDissipation(const viscousDissipation&) = delete;

        //- No copy assignment
        void operator=(const viscousDissipation&) = delete;


    //- Destructor
    virtual ~viscousDissipation() = default;


    // Member Functions

        //- Add explicit contribution to compressible energy equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );


        //- Read source dictionary - no-op
        virtual bool read(const dictionary& dict)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

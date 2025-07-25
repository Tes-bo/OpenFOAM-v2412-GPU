/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::buoyancyForce

Group
    grpFvOptionsSources

Description
    Applies buoyancy force sources (i.e. \c Sf) to the momentum
    equation for incompressible and compressible flows.

    \f[
        S_f = \rho \vec{g}
    \f]

    where
    \vartable
      S_f     | Buoyancy force source term
      \rho    | Fluid density
      \vec{g} | Gravitational acceleration field
    \endvartable

    Sources applied to:
    \verbatim
      U         | Velocity                                   [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity                                   [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    buoyancyForce1
    {
        // Mandatory entries (unmodifiable)
        type             buoyancyForce;
        fields           (U);

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: buoyancyEnergy         | word   | yes   | -
      fields    | Name of operand vector field      | word   | yes   | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

Note
  - The entry \c fields accepts only a single vector field name.

See also
  - Foam::fv::buoyancyEnergy
  - Foam::fv::buoyancyTurbSource

SourceFiles
    buoyancyForce.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyForce_H
#define buoyancyForce_H

#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class buoyancyForce Declaration
\*---------------------------------------------------------------------------*/

class buoyancyForce
:
    public fv::option
{
    // Private Data

        const uniformDimensionedVectorField g_;


public:

    //- Runtime type information
    TypeName("buoyancyForce");


    // Constructors

        //- Construct from explicit source name and mesh
        buoyancyForce
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        buoyancyForce(const buoyancyForce&) = delete;

        //- No copy assignment
        void operator=(const buoyancyForce&) = delete;


    //- Destructor
    virtual ~buoyancyForce() = default;


    // Member Functions

        //- Add explicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );


        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

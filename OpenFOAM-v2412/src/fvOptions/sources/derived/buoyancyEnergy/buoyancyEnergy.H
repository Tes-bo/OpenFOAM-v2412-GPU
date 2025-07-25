/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::fv::buoyancyEnergy

Group
    grpFvOptionsSources

Description
    Applies buoyancy energy sources (i.e. \c Sb)
    to the energy equation for compressible flows.

    \f[
        S_b = \rho (\vec{U} \cdot \vec{g})
    \f]

    where
    \vartable
      S_b     | Buoyancy energy source term
      \rho    | Fluid density
      \vec{U} | Velocity field
      \vec{g} | Gravitational acceleration field
    \endvartable

    Sources applied to either of the below, if exist:
    \verbatim
      e         | Internal energy                            [m2/s2]
      h         | Enthalphy                                  [m2/s2]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity                                   [m/s]
      rho       | Fluid density                              [kg/m3]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    buoyancyEnergy1
    {
        // Mandatory entries (unmodifiable)
        type             buoyancyEnergy;
        fields           (h);

        // Optional entries (runtime modifiable)
        U                U;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: buoyancyEnergy         | word   | yes   | -
      fields    | Name of operand energy field      | word   | yes   | -
      U         | Name of operand velocity field    | word   | no    | U
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

Note
  - The entry \c fields accepts only a single field name.
  - The sources are applied onto the entire computational domain.

See also
  - Foam::fv::buoyancyForce
  - Foam::fv::buoyancyTurbSource

SourceFiles
    buoyancyEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancyEnergy_H
#define buoyancyEnergy_H

#include "fvOption.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class buoyancyEnergy Declaration
\*---------------------------------------------------------------------------*/

class buoyancyEnergy
:
    public fv::option
{
    // Private Data

        //- Name of operand velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("buoyancyEnergy");


    // Constructors

        //- Construct from explicit source name and mesh
        buoyancyEnergy
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        buoyancyEnergy(const buoyancyEnergy&) = delete;

        //- No copy assignment
        void operator=(const buoyancyEnergy&) = delete;


    //- Destructor
    virtual ~buoyancyEnergy() = default;


    // Member Functions

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
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

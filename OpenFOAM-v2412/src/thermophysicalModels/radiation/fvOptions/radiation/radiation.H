/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::fv::radiation

Description
    Applies radiation sources (i.e. \c Sh)
    to the energy equation for compressible flows.

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    radiation1
    {
        // Mandatory entries (unmodifiable)
        type             radiation;
        fields           (h);

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: radiation              | word   | yes   | -
      fields    | Name of operand energy field      | word   | yes   | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink

SourceFiles
    radiation.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_H
#define radiation_H

#include "fvOption.H"
#include "uniformDimensionedFields.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class radiation Declaration
\*---------------------------------------------------------------------------*/

class radiation
:
    public fv::option
{
    // Private Data

        //- The radiation model pointer
        autoPtr<Foam::radiation::radiationModel> radiation_;


    // Private Member Functions

        //- No copy construct
        radiation(const radiation&) = delete;

        //- No copy assignment
        void operator=(const radiation&) = delete;


public:

    //- Runtime type information
    TypeName("radiation");


    // Constructors

        //- Construct from explicit source name and mesh
        radiation
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Evaluate

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

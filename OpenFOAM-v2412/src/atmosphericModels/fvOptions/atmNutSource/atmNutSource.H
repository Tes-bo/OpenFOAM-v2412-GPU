/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
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
    Foam::fv::atmNutSource

Group
    grpFvOptionsSources

Description
    Adds/subtracts a given artificial turbulent viscosity field to/from \c nut
    for atmospheric boundary layer modelling.

    Corrections applied to:
    \verbatim
      nut         | Turbulent viscosity                  [m2/s]
    \endverbatim

    Required fields:
    \verbatim
      nut         | Turbulent viscosity                  [m2/s]
      artNut      | Artificial turbulent viscosity       [m2/s]
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    atmNutSource1
    {
        // Mandatory entries (unmodifiable)
        type             atmNutSource;

        atmNutSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode    all;

            // Optional entries (unmodifiable)
            nut              artNut;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property      | Description                  | Type   | Req'd  | Dflt
      type          | Type name: atmNutSource      | word   | yes    | -
      nut           | Field name of artificial nut | word   | no     | artNut
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

SourceFiles
    atmNutSource.C

\*---------------------------------------------------------------------------*/

#ifndef fv_atmNutSource_H
#define fv_atmNutSource_H

#include "cellSetOption.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class atmNutSource Declaration
\*---------------------------------------------------------------------------*/

class atmNutSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Name of the artificial turbulent viscosity field
        const word artNutName_;


    // Fields

        //- Artificial turbulent viscosity field [m2/s]
        volScalarField artNut_;


public:

    //- Runtime type information
    TypeName("atmNutSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmNutSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmNutSource(const atmNutSource&) = delete;

        //- No copy assignment
        void operator=(const atmNutSource&) = delete;


    // Member Functions

        //- Correct the turbulent viscosity
        virtual void correct(volScalarField& field);

        //- Read source dictionary (effectively no-op)
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

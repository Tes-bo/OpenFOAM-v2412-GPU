/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::fv::variableHeatTransfer

Group
    grpFvOptionsSources

Description
    Applies a variable heat transfer model depending
    on local values for inter-region heat exchanges.

    The area of contact between regions (area) must be provided.
    The Nu number is calculated as:

    \f[
        Nu = a \, \textrm{Re}^b \, \textrm{Pr}^c
    \f]

    and the heat transfer coefficient as:

    \f[
        htc = \frac{\textrm{Nu} \, K}{ds}
    \f]

    where
    \vartable
      htc     | Heat transfer coefficient
      Nu      | Nusselt number
      Re      | Reynolds number
      Pr      | Fluid Prandtl number
      K       | Heat conduction
      ds      | Strut diameter
      a, b, c | Model constants
    \endvartable

    Required fields:
    \verbatim
      AoV       | Area per unit volume of heat exchanger    [1/m]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    variableHeatTransfer1
    {
        // Mandatory entries (unmodifiable)
        type             variableHeatTransfer;

        // Optional entries (runtime modifiable)
        UNbr             U;
        a                0.0;
        b                0.0;
        c                0.0;
        ds               0.0;
        Pr               0.0;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                     | Type  | Reqd | Dflt
      type      | Type name: variableHeatTransfer | word  | yes  | -
      UNbr      | Name of operand neighbour velocity field | word | no | U
      a         | Model constant                  | scalar | no  | 0
      b         | Model constant                  | scalar | no  | 0
      c         | Model constant                  | scalar | no  | 0
      ds        | Struct diameter                 | scalar | no  | 0
      Pr        | Fluid Prandtl number            | scalar | no  | 0
      \endtable

    The inherited entries are elaborated in:
     - \link interRegionHeatTransfer.H \endlink

See also
  - Foam::fv::interRegionHeatTransferModel
  - Foam::fv::tabulatedHeatTransfer
  - Foam::fv::tabulatedNTUHeatTransfer
  - Foam::fv::constantHeatTransfer

SourceFiles
    variableHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef variableHeatTransfer_H
#define variableHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class variableHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class variableHeatTransfer
:
    public interRegionHeatTransferModel
{

    // Private Data

        //- Name of operand neighbour velocity field
        word UNbrName_;

        //- Model constants
        scalar a_;
        scalar b_;
        scalar c_;

        //- Strut diameter
        scalar ds_;

        //- Fluid Prandtl number
        scalar Pr_;

        //- Area per unit volume of heat exchanger
        autoPtr<volScalarField> AoV_;


public:

    //- Runtime type information
    TypeName("variableHeatTransfer");


    // Constructors

        //- Construct from dictionary
        variableHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        variableHeatTransfer(const variableHeatTransfer&) = delete;

        //- No copy assignment
        void operator=(const variableHeatTransfer&) = delete;


    //- Destructor
    virtual ~variableHeatTransfer() = default;


    // Public Functions

        //- Calculate the heat transfer coefficient
        virtual void calculateHtc();

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

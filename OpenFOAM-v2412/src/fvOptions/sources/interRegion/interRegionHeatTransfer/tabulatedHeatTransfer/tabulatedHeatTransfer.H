/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fv::tabulatedHeatTransfer

Group
    grpFvOptionsSources

Description
    Applies a tabulated heat transfer model for inter-region heat exchanges.

    The 2D table returns the heat transfer coefficient
    by querying the local and neighbour region velocities.

    Required fields:
    \verbatim
      AoV       | Area per unit volume of heat exchanger    [1/m]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    tabulatedHeatTransfer1
    {
        // Mandatory entries (unmodifiable)
        type            tabulatedHeatTransfer;
        outOfBounds     clamp;
        file            "table";

        // Optional entries (unmodifiable)
        U               <Uname>;
        UNbr            <UnbrName>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                      | Type  | Reqd | Dflt
      type      | Type name: tabulatedHeatTransfer | word  | yes  | -
      file  | Heat transfer coefficient table  | interpolate2DTable | yes | -
      U         | Name of operand velocity field   | word  | no   | U
      UNbr      | Name of operand neighbour velocity field | word | no | U
    \endtable

    The inherited entries are elaborated in:
     - \link interRegionHeatTransfer.H \endlink
     - \link interpolate2DTable.H \endlink

See also
  - Foam::fv::interRegionHeatTransferModel
  - Foam::fv::constantHeatTransfer
  - Foam::fv::tabulatedNTUHeatTransfer
  - Foam::fv::variableHeatTransfer

SourceFiles
    tabulatedHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef tabulatedHeatTransfer_H
#define tabulatedHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class tabulatedHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class tabulatedHeatTransfer
:
    public interRegionHeatTransferModel
{

    // Private Data

        //- Name of operand velocity field
        word UName_;

        //- Name of neighbour velocity field
        word UNbrName_;

        //- 2D look up table
        autoPtr<interpolation2DTable<scalar>> hTable_;

        //- Area per unit volume of heat exchanger
        autoPtr<volScalarField> AoV_;

        //- Heat transfer coefficient table
        const interpolation2DTable<scalar>& hTable();

        //- Field of area divided by volume
        const volScalarField& AoV();

        //- Start time name
        const word startTimeName_;


public:

    //- Runtime type information
    TypeName("tabulatedHeatTransfer");


    // Constructors

        //- Construct from dictionary
        tabulatedHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        tabulatedHeatTransfer(const tabulatedHeatTransfer&) = delete;

        //- No copy assignment
        void operator=(const tabulatedHeatTransfer&) = delete;


    //- Destructor
    virtual ~tabulatedHeatTransfer() = default;


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

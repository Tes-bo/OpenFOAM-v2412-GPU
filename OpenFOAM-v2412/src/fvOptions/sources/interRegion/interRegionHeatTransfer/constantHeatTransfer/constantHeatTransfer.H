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
    Foam::fv::constantHeatTransfer

Group
    grpFvOptionsSources

Description
    Applies a constant heat transfer model for inter-region heat exchanges.

    Required fields:
    \verbatim
      htcConst  | Constant heat transfer coefficient        [W/m2/K]
      AoV       | Area per unit volume of heat exchanger    [1/m]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    constantHeatTransfer1
    {
        // Mandatory entries (unmodifiable)
        type                 constantHeatTransfer;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                     | Type  | Reqd | Dflt
      type      | Type name: constantHeatTransfer | word  | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link interRegionHeatTransfer.H \endlink

See also
  - Foam::fv::interRegionHeatTransferModel
  - Foam::fv::tabulatedHeatTransfer
  - Foam::fv::tabulatedNTUHeatTransfer
  - Foam::fv::variableHeatTransfer

SourceFiles
    constantHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef constantHeatTransfer_H
#define constantHeatTransfer_H

#include "interRegionHeatTransferModel.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class constantHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class constantHeatTransfer
:
    public interRegionHeatTransferModel
{
    // Private Data

        //- Constant heat transfer coefficient [W/m2/K]
        autoPtr<volScalarField> htcConst_;

        //- Area per unit volume of heat exchanger [1/m]
        autoPtr<volScalarField> AoV_;


public:

    //- Runtime type information
    TypeName("constantHeatTransfer");


    // Constructors

        //- Construct from components
        constantHeatTransfer
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        constantHeatTransfer(const constantHeatTransfer&) = delete;

        //- No copy assignment
        void operator=(const constantHeatTransfer&) = delete;


    //- Destructor
    virtual ~constantHeatTransfer() = default;


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

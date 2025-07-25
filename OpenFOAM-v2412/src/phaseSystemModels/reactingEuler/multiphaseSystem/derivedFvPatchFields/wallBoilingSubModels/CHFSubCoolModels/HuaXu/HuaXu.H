/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd
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
    Foam::wallBoilingModels::CHFModels::HuaXu

Description
    A critical heat flux (CHF) sub-cooling correlation model
    based on Hua-Xu (2000) for boiling flows.

    Reference:
    \verbatim
        Hua, T. C., & Xu, J. J. (2000).
        Quenching boiling in subcooled liquid nitrogen
        for solidification of aqueous materials.
        Materials Science and Engineering: A, 292(2), 169-172.
        DOI:10.1016/S0921-5093(00)01004-2
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    CHFSubCoolModel
    {
        // Mandatory entries
        type            HuaXu;

        // Optional entries
        Kburn           <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: HuaXu                    | word | yes  | -
      Kburn     | Burn out factor                   | scalar | no   | 1.5
    \endtable

SourceFiles
    HuaXu.C

\*---------------------------------------------------------------------------*/

#ifndef HuaXu_H
#define HuaXu_H

#include "CHFSubCoolModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace CHFModels
{

/*---------------------------------------------------------------------------*\
                         Class HuaXu Declaration
\*---------------------------------------------------------------------------*/

class HuaXu
:
    public CHFSubCoolModel
{
    // Private Data

        //- Burn out factor
        scalar Kburn_;


    // Private Member Functions

        //- No copy construct
        HuaXu(const HuaXu&) = delete;

        //- No copy assignment
        void operator=(const HuaXu&) = delete;


public:

    //- Runtime type information
    TypeName("HuaXu");


    // Constructors

        //- Construct from a dictionary
        HuaXu(const dictionary& dict);


    //- Destructor
    virtual ~HuaXu() = default;


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> CHFSubCool
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace CHFModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

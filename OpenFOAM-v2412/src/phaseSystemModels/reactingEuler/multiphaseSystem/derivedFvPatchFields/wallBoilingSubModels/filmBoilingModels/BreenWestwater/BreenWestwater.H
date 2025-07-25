/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd
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
    Foam::wallBoilingModels::filmBoilingModels::BreenWestwater

Description
    Boiling film correlation.
    A correlation for boiling film modelling
    based on Breen & Westwater (1965) for boiling flows.

    References:
    \verbatim
        Breen B.P. & Westwater J. W. (1965)
        Effect of diameter of horizontal
        tubes on film boiling heat tranfer.
        Chem. Eng. Progr. 58. No 7.
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    filmBoilingModel
    {
        // Mandatory entries
        type            BreenWestwater;

        // Optional entries
        Cn              <scalar>;
        an              <scalar>;
        bn              <scalar>;
        n               <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                       | Type | Reqd | Deflt
      type       | Type name: BreenWestwater         | word | yes  | -
      Cn         | Model coefficient               | scalar | no   | 0.37
    \endtable

SourceFiles
    BreenWestwater.C

\*---------------------------------------------------------------------------*/

#ifndef BreenWestwater_H
#define BreenWestwater_H

#include "filmBoilingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace filmBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class BreenWestwater Declaration
\*---------------------------------------------------------------------------*/

class BreenWestwater
:
    public filmBoilingModel
{
    // Private Data

        //- Model coefficient
        scalar Cn_;


    // Private Member Functions

        //- No copy construct
        BreenWestwater(const BreenWestwater&) = delete;

        //- No copy assignment
        void operator=(const BreenWestwater&) = delete;


public:

    //- Runtime type information
    TypeName("BreenWestwater");


    // Constructors

        //- Construct from a dictionary
        BreenWestwater(const dictionary& dict);


    //- Destructor
    virtual ~BreenWestwater() = default;


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> htcFilmBoil
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace filmBoilingModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

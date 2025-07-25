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
    Foam::wallBoilingModels::filmBoilingModels::Bromley

Description
    A correlation for boiling film modelling
    based on Bromley (1950) for boiling flows.

    Reference:
    \verbatim
        Bromley, L. A. (1950).
        Heat transfer in stable film boiling.
        Chemical Engineering Progress, 46, 221-227.
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    filmBoilingModel
    {
        // Mandatory entries
        type            Bromley;
        L               <scalar>;

        // Optional entries
        Cn              <scalar>;
        emissivity      <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: Bromley                  | word | yes  | -
      L         | Characteristic length scale       | scalar | yes  | -
      Cn  | Coefficient for nucleation site density | scalar | no   | 0.62
      emissivity | Wall emissivity                  | scalar | no   | 1.0
    \endtable

SourceFiles
    Bromley.C

\*---------------------------------------------------------------------------*/

#ifndef Bromley_H
#define Bromley_H

#include "filmBoilingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace filmBoilingModels
{

/*---------------------------------------------------------------------------*\
                         Class Bromley Declaration
\*---------------------------------------------------------------------------*/

class Bromley
:
    public filmBoilingModel
{
    // Private Data

        //- Coefficient for nucleation site density
        scalar Cn_;

        //- Wall emissivity
        scalar emissivity_;

        //- Characteristic length scale
        scalar L_;


    // Private Member Functions

        //- No copy construct
        Bromley(const Bromley&) = delete;

        //- No copy assignment
        void operator=(const Bromley&) = delete;


public:

    //- Runtime type information
    TypeName("Bromley");

    // Constructors

        //- Construct from a dictionary
        Bromley(const dictionary& dict);


    //- Destructor
    virtual ~Bromley() = default;


    // Member Functions

        //- Calculate and return the film boiling correlation
        virtual tmp<scalarField> htcFilmBoil
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

} // End namespace filmBoilingModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

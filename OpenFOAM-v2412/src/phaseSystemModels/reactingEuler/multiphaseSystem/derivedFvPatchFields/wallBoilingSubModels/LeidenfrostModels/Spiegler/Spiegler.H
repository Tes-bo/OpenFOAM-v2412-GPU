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
    Foam::wallBoilingModels::LeidenfrostModels::Spiegler

Description
    A model for Leidenfrost effects based on
    Spiegler et al. (1963) for boiling flows.

    Reference:
    \verbatim
        Spiegler, P., Hopenfeld, J., Silberberg, M.,
        Bumpus Jr, C. F., & Norman, A. (1963).
        Onset of stable film boiling and the foam limit.
        International Journal of Heat and Mass Transfer, 6(11), 987-989.
        DOI:10.1016/0017-9310(63)90053-X
    \endverbatim

Usage
    Example of the model specification:
    \verbatim
    LeidenfrostModel
    {
        // Mandatory entries
        type            Spiegler;

        // Optional entries
        Tcrit            <scalar>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type | Reqd | Deflt
      type      | Type name: Spiegler                 | word | yes  | -
      Tcrit     | Critical temperature [K]          | scalar | no   | 374.0
    \endtable

SourceFiles
    Spiegler.C

\*---------------------------------------------------------------------------*/

#ifndef Spiegler_H
#define Spiegler_H

#include "LeidenfrostModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace LeidenfrostModels
{

/*---------------------------------------------------------------------------*\
                         Class Spiegler Declaration
\*---------------------------------------------------------------------------*/

class Spiegler
:
    public LeidenfrostModel
{
    // Private Data

        //- Critical temperature
        scalar Tcrit_;


    // Private Member Functions

        //- No copy construct
        Spiegler(const Spiegler&) = delete;

        //- No copy assignment
        void operator=(const Spiegler&) = delete;


public:

    //- Runtime type information
    TypeName("Spiegler");


    // Constructors

        //- Construct from a dictionary
        Spiegler(const dictionary& dict);


    //- Destructor
    virtual ~Spiegler() = default;


    // Member Functions

        //- Calculate and return the Leidenfrost temperature
        virtual tmp<scalarField> TLeid
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

} // End namespace LeidenfrostModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

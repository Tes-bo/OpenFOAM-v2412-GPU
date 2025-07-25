/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::solidArrheniusReactionRate

Description
    Arrhenius reaction rate for solids

SourceFiles
    solidArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_solidArrheniusReactionRate_H
#define Foam_solidArrheniusReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class solidArrheniusReactionRate;

Ostream& operator<<(Ostream&, const solidArrheniusReactionRate&);


/*---------------------------------------------------------------------------*\
                  Class solidArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class solidArrheniusReactionRate
{
    // Private Data

        //- Pre-exponential factor
        scalar A_;

        //- Activation temperature
        scalar Ta_;

        //- Critical temperature
        scalar Tcrit_;


public:

    // Constructors

        //- Construct from components
        inline solidArrheniusReactionRate
        (
            const scalar A,
            const scalar Ta,
            const scalar Tcrit
        );


        //- Construct from dictionary
        inline solidArrheniusReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    //- Destructor
    virtual ~solidArrheniusReactionRate() = default;


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "Arrhenius";
        }

        //- Return reaction rate constant
        //  frequency of collisions resulting in a reaction
        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;


        //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const solidArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::LESdelta

Description
    Abstract base class for LES deltas

SourceFiles
    LESdelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESdelta_H
#define LESdelta_H

#include "turbulenceModel.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class LESdelta Declaration
\*---------------------------------------------------------------------------*/

class LESdelta
{
protected:

    // Protected Data

        const turbulenceModel& turbulenceModel_;

        volScalarField delta_;


    // Protected Member Functions

        //- No copy construct
        LESdelta(const LESdelta&) = delete;

        //- No copy assignment
        void operator=(const LESdelta&) = delete;


public:

    //- Runtime type information
    TypeName("LESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESdelta,
            dictionary,
            (
                const word& name,
                const turbulenceModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        LESdelta
        (
            const word& name,
            const turbulenceModel& turbulence
        );


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary& dict,
            const word& lookupName = "delta"
        );

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTableType& additionalConstructors,
            const word& lookupName = "delta"
        );


    //- Destructor
    virtual ~LESdelta() = default;


    // Member Functions

        //- Return turbulenceModel reference
        const turbulenceModel& turbulence() const
        {
            return turbulenceModel_;
        }

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        virtual operator const volScalarField&() const
        {
            return delta_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

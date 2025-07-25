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

Namespace
    Foam::solidBodyMotionFunctions

Description
    Namespace for solid-body motions


Class
    Foam::solidBodyMotionFunction

Description
    Base class for defining solid-body motions

SourceFiles
    solidBodyMotionFunction.C
    solidBodyMotionFunctionNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_solidBodyMotionFunction_H
#define Foam_solidBodyMotionFunction_H

#include "Time.H"
#include "dictionary.H"
#include "septernion.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solidBodyMotionFunction Declaration
\*---------------------------------------------------------------------------*/

class solidBodyMotionFunction
{
protected:

    // Protected Data

        dictionary SBMFCoeffs_;

        const Time& time_;


    // Protected Member Functions

        //- No copy construct
        solidBodyMotionFunction(const solidBodyMotionFunction&) = delete;

        //- No copy assignment
        void operator=(const solidBodyMotionFunction&) = delete;


public:

    //- Runtime type information
    TypeName("solidBodyMotionFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidBodyMotionFunction,
            dictionary,
            (const dictionary& SBMFCoeffs, const Time& runTime),
            (SBMFCoeffs, runTime)
        );


    // Constructors

        //- Construct from the coefficients dictionary and Time
        solidBodyMotionFunction
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const = 0;


    // Selectors

        //- Construct and dispatch motionType with dictionary and Time.
        //  Returns nullptr if motionType is empty
        static autoPtr<solidBodyMotionFunction> New
        (
            const word& motionType,
            const dictionary& dict,
            const Time& runTime
        );

        //- Select "solidBodyMotionFunction" type from dictionary
        //- and create with Time.
        static autoPtr<solidBodyMotionFunction> New
        (
            const dictionary& dict,
            const Time& runTime
        );

        //- Select "solidBodyMotionFunction" type (if present) from dictionary
        //- and create with Time.
        static autoPtr<solidBodyMotionFunction> NewIfPresent
        (
            const dictionary& dict,
            const Time& runTime
        );


    //- Destructor
    virtual ~solidBodyMotionFunction() = default;


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const = 0;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs) = 0;

        //- Write in dictionary format
        virtual void writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

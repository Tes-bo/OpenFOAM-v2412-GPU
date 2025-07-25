/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::bisection

Description
    Reduces step by a given ratio

SourceFiles
    bisection.C

\*---------------------------------------------------------------------------*/

#ifndef bisection_H
#define bisection_H

#include "stepUpdate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class bisection Declaration
\*---------------------------------------------------------------------------*/

class bisection
:
    public stepUpdate
{
protected:

    // Protected data

        scalar ratio_;


private:

    // Private Member Functions

        //- No copy construct
        bisection(const bisection&) = delete;

        //- No copy assignment
        void operator=(const bisection&) = delete;


public:

    //- Runtime type information
    TypeName("bisection");


    // Constructors

        //- Construct from components
        bisection(const dictionary& dict);


    //- Destructor
    virtual ~bisection() = default;


    // Member Functions

       // Update step in given direction
       virtual void updateStep(scalar& step);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

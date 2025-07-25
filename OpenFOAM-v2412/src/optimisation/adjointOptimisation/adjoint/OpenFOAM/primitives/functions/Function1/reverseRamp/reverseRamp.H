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
    Foam::Function1Types::reverseRamp

Description
    Reverse ramp function starting from 1 and decreasing stepRamply to a min. 
    user-specified value from \c start over the \c duration and remaining 
    constant thereafter.

See also 
    Foam::Function1Types::ramp

SourceFiles
    reverseRamp.C

\*---------------------------------------------------------------------------*/

#ifndef reverseRamp_H
#define reverseRamp_H

#include "ramp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                        Class reverseRamp Declaration
\*---------------------------------------------------------------------------*/

class reverseRamp
:
    public ramp
{
protected:

    // Protected Data Members 

        //- The min value 
        scalar minValue_;

        //- Interval for reducing the base value
        scalar interval_;

        //- Steps to reach the min. value
        scalar steps_;


private:

    // Private Data Members 

        //- No copy assignment
        void operator=(const reverseRamp&) = delete;


public:

    //- Runtime type information
    TypeName("reverseRamp");


    // Constructors

        //- Construct from entry and dictionary 
        reverseRamp 
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );


    //- Destructor
    virtual ~reverseRamp() = default;


    // Member Functions

        //- Return value for time t
        inline virtual scalar value(const scalar t) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reverseRampI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

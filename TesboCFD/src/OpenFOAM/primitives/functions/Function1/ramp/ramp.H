/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::Function1Types::ramp

Description
    Ramp function base class for the set of scalar functions starting from 0 and
    increasing monotonically to 1 from \c start over the \c duration and
    remaining at 1 thereafter.

    Usage:
    \verbatim
        <entryName> <rampFunction>;
        <entryName>Coeffs
        {
            start     10;
            duration  20;
        }
    \endverbatim
    or
    \verbatim
        <entryName>
        {
            type      <rampFunction>;
            start     10;
            duration  20;
        }
    \endverbatim

    Where:
    \table
        Property | Description  | Required | Default value
        start    | Start time   | no       | 0
        duration | Duration     | yes      |
    \endtable

See also
    Foam::Function1

SourceFiles
    ramp.C

\*---------------------------------------------------------------------------*/

#ifndef Function1Types_ramp_H
#define Function1Types_ramp_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           Class ramp Declaration
\*---------------------------------------------------------------------------*/

class ramp
:
    public Function1<scalar>
{
protected:

    // Protected Data

        //- Start-time of the ramp function
        scalar start_;

        //- Duration of the ramp function
        scalar duration_;

        //- Simple linear ramp function
        //- that forms the basis of many more complex ramp functions
        inline scalar linearRamp(const scalar t) const
        {
            return clamp((t - start_)/duration_, zero_one{});
        }


private:

    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);

        //- No copy assignment
        void operator=(const ramp&) = delete;


public:

    // Constructors

        //- Construct from entry name, dictionary and optional registry
        ramp
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );


    //- Destructor
    virtual ~ramp() = default;


    // Member Functions

        //- Convert time
        virtual void userTimeToTime(const Time& t);

        //- Return value for time t
        virtual scalar value(const scalar t) const = 0;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        virtual void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

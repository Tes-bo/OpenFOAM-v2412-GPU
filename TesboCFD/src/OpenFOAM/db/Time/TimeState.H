/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::TimeState

Description
    The time value with time-stepping information, user-defined remapping, etc.

SourceFiles
    TimeState.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_TimeState_H
#define Foam_TimeState_H

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TimeState Declaration
\*---------------------------------------------------------------------------*/

class TimeState
:
    public dimensionedScalar
{
    // Private Member Functions

        // //- Dimensions for 1/time
        // inline static dimensionSet dimInvTime()
        // {
        //     return dimensionSet(0, 0, -1, 0, 0, 0, 0);
        // }


protected:

    // Protected Data

        label timeIndex_;
        label writeTimeIndex_;

        scalar deltaT_;
        scalar deltaT0_;
        scalar deltaTSave_;

        bool deltaTchanged_;
        bool writeTime_;


    // Protected Member Functions

        //- Reset some of TimeState (name, value, index)
        void resetTimeState
        (
            const word& newName,
            const scalar newValue,
            const label newIndex
        )
        {
            dimensionedScalar::name() = newName;
            dimensionedScalar::value() = newValue;
            timeIndex_ = newIndex;
        }


public:

    // Constructors

        //- Construct a zero time state, using the current time formatting
        TimeState();


    //- Destructor
    virtual ~TimeState() = default;


    // Member Functions

    // Conversion

        //- Convert the user-time (e.g. CA deg) to real-time (s).
        virtual scalar userTimeToTime(const scalar theta) const;

        //- Convert the real-time (s) into user-time (e.g. CA deg)
        virtual scalar timeToUserTime(const scalar t) const;


    // Access

        //- Return the current user-time value.
        //- (ie, after applying any timeToUserTime() conversion)
        inline scalar timeOutputValue() const;

        //- Return the current time name
        inline const word& timeName() const noexcept;

        //- Return the current time index
        inline label timeIndex() const noexcept;

        //- Return time step value
        inline scalar deltaTValue() const noexcept;

        //- Return old time step value
        inline scalar deltaT0Value() const noexcept;

        //- Return time step
        inline dimensionedScalar deltaT() const;

        //- Return old time step
        inline dimensionedScalar deltaT0() const;

        //- True if this is a write interval
        inline bool writeTime() const noexcept;


    // Housekeeping

        //- Deprecated(2016-05) return true if this is a write time.
        //  \deprecated(2016-05) - use writeTime() method
        bool outputTime() const noexcept { return this->writeTime(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TimeStateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

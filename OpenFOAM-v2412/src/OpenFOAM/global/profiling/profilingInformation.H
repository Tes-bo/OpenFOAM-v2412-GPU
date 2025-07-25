/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2016 Bernhard Gschaider
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::profilingInformation

Description
    Code profiling information in terms of time spent, number of calls etc.

SourceFiles
    profilingInformation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_profilingInformation_H
#define Foam_profilingInformation_H

#include "labelFwd.H"
#include "scalarFwd.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class profiling;

/*---------------------------------------------------------------------------*\
                    Class profilingInformation Declaration
\*---------------------------------------------------------------------------*/

class profilingInformation
{
    // Private Data

        //- Unique id to identify it
        const label id_;

        //- What this timer does
        const std::string description_;

        //- Pointer to the parent object (or self for top-level)
        profilingInformation* parent_;

        //- Number of times this was called
        long calls_;

        //- Total time spent
        scalar totalTime_;

        //- Time spent in children
        scalar childTime_;

        //- Max memory usage on call.
        //  Only valid when the calling profiling has memInfo active.
        mutable int maxMem_;

        //- Is this information active or passive (ie, on the stack)?
        mutable bool active_;


protected:

    // Friendship

        //- Allow creation of master-element, setActive etc.
        friend class profiling;


    // Constructors

        //- Default construct - only the master-element
        profilingInformation();

        //- No copy construct
        profilingInformation(const profilingInformation&) = delete;

        //- No copy assignment
        void operator=(const profilingInformation&) = delete;


    // Protected Member Functions

        //- Mark as being active or passive
        void setActive(bool state) const;

public:

    // Constructors

        //- Construct from components
        profilingInformation
        (
            profilingInformation* parent,
            const std::string& descr,
            const label id
        );


    //- Destructor
    ~profilingInformation() = default;


    // Member Functions

    // Access

        label id() const noexcept { return id_; }

        const std::string& description() const noexcept { return description_; }

        profilingInformation& parent() const noexcept { return *parent_; }

        label calls() const noexcept { return calls_; }

        scalar totalTime() const noexcept { return totalTime_; }

        scalar childTime() const noexcept { return childTime_; }

        int maxMem() const noexcept { return maxMem_; }

        bool active() const noexcept { return active_; }


    // Edit

        //- Update it with a new timing information
        void update(const scalar elapsedTime);


    // Write

        //- Write the profiling times, optionally with additional values
        //  Uses dictionary format.
        Ostream& write
        (
            Ostream& os,
            const bool offset = false,
            const scalar elapsedTime = 0,
            const scalar childTime = 0
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators

Ostream& operator<<(Ostream& os, const profilingInformation& info);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

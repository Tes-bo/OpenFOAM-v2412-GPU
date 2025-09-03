/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::memInfo

Description
    Memory usage information for the current process, and the system memory
    that is free.

Note
    Uses the information from /proc/PID/status and from /proc/meminfo

SourceFiles
    memInfo.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_memInfo_H
#define Foam_memInfo_H

#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class memInfo;
class word;
class Ostream;

Ostream& operator<<(Ostream& os, const memInfo& m);

/*---------------------------------------------------------------------------*\
                           Class memInfo Declaration
\*---------------------------------------------------------------------------*/

class memInfo
{
    // Private Data

        //- Peak memory used by the process (VmPeak in /proc/PID/status)
        int64_t peak_;

        //- Memory used by the process (VmSize in /proc/PID/status)
        int64_t size_;

        //- Resident set size of the process (VmRSS in /proc/PID/status)
        int64_t rss_;

        //- System memory free (MemFree in /proc/meminfo)
        int64_t free_;


    // Private Member Functions

        //- Populate entries
        void populate();


public:

    // Constructors

        //- Construct and populate with values
        memInfo();


    // Member Functions

    // Access

        //- True if the memory information appears valid
        bool good() const noexcept;

        //- Peak memory at last update - (VmPeak in /proc/PID/status)
        int64_t peak() const noexcept { return peak_; }

        //- Memory size at last update - (VmSize in /proc/PID/status)
        int64_t size() const noexcept { return size_; }

        //- Resident set size at last update - (VmRSS in /proc/PID/status)
        int64_t rss() const noexcept { return rss_; }

        //- System memory free (MemFree in /proc/meminfo)
        int64_t free() const noexcept { return free_; }


    // Edit

        //- Reset to zero
        void clear() noexcept;

        //- Update according to /proc/PID/status and /proc/memory contents
        const memInfo& update();


    // IOstream Operators

        //- Write peak/size/rss to stream
        friend Ostream& operator<<(Ostream& os, const memInfo& m);


    // Write

        //- Write mem-info as dictionary entries
        void writeEntries(Ostream& os) const;

        //- Write mem-info as dictionary
        void writeEntry(const word& keyword, Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

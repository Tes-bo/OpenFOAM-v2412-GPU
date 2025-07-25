/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::threadedCollatedOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    threadedCollatedOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_threadedCollatedOFstream_H
#define Foam_threadedCollatedOFstream_H

#include "dictionary.H"
#include "StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class OFstreamCollator;

/*---------------------------------------------------------------------------*\
                  Class threadedCollatedOFstream Declaration
\*---------------------------------------------------------------------------*/

class threadedCollatedOFstream
:
    public OStringStream
{
    // Private Data

        //- The backend writer
        OFstreamCollator& writer_;

        //- The backend file name
        const fileName pathName_;

        //- Atomic file creation (ignored with append)
        const IOstreamOption::atomicType atomic_;

        //- Output file compression
        const IOstreamOption::compressionType compression_;

        //- Use threading
        const bool useThread_;

        //- Additional FoamFile entries for decomposed data
        dictionary headerEntries_;


public:

    // Constructors

        //- Construct and set stream status
        threadedCollatedOFstream
        (
            OFstreamCollator& writer,
            IOstreamOption::atomicType atomic,
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            const bool useThread = true
        );

        //- Construct and set stream status
        threadedCollatedOFstream
        (
            OFstreamCollator& writer,
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            const bool useThread = true
        );


    //- Destructor
    ~threadedCollatedOFstream();


    // Member Functions

        //- Define the header entries for the data block(s)
        void setHeaderEntries(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

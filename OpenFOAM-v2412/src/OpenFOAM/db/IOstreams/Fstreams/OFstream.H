/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::OFstream

Description
    Output to file stream as an OSstream, normally using \c std::ofstream
    for the actual output.

Note
    The atomic output works by creating an intermediate temporary file,
    which is renamed as an atomic operation when closing. It is not
    possible, or particularly desirable, to have an atomic in combination
    with append behaviour. If both are specified, append has priority.

Note
    An output file can be opened in two different \c append modes, both of
    which preserve existing files:
    -# A common append mode is APPEND_APP, which corresponds to the
       \c std::ios_base::app flag.
       A seek-to-end is performed at \em every write.
       It is thus not possible to use any manual seeks to overwrite parts
       of the file.
    -# The other append mode is APPEND_ATE, which roughly corresponds to the
       \c std::ios_base::ate flag behaviour.
       A seek-to-end is performed immediately after opening,
       but not subsequently.
       Manual seeks can be used to overwrite parts of the file.
    .

SourceFiles
    OFstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OFstream_H
#define Foam_OFstream_H

#include "OSstream.H"
#include "className.H"
#include "fstreamPointer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OFstream Declaration
\*---------------------------------------------------------------------------*/

class OFstream
:
    private Foam::ofstreamPointer,
    public OSstream
{
public:

    //- Declare type-name (with debug switch)
    ClassName("OFstream");


    // Constructors

        //- Construct a null output file stream that behaves like \c /dev/null
        explicit OFstream(std::nullptr_t);

        //- Construct with specified atomic behaviour
        //- from pathname, stream option, optional append (see note).
        OFstream
        (
            IOstreamOption::atomicType atomic,
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            IOstreamOption::appendType append = IOstreamOption::NO_APPEND
        );

        //- Construct from pathname and other specifications.
        //  See note on append mode.
        explicit OFstream
        (
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            IOstreamOption::appendType append = IOstreamOption::NO_APPEND
        )
        :
            OFstream(IOstreamOption::NON_ATOMIC, pathname, streamOpt, append)
        {}

        //- Construct from pathname, format (uncompressed),
        //- optional append (see note),
        //- atomic behaviour as per system default
        OFstream
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED,
            IOstreamOption::appendType append = IOstreamOption::NO_APPEND
        )
        :
            OFstream(pathname, IOstreamOption(fmt, cmp), append)
        {}

        //- Construct with specified atomic behaviour
        //- from pathname, format (uncompressed),
        //- optional append (see note).
        OFstream
        (
            IOstreamOption::atomicType atomic,
            const fileName& pathname,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED,
            IOstreamOption::appendType append = IOstreamOption::NO_APPEND
        )
        :
            OFstream(atomic, pathname, IOstreamOption(fmt, cmp), append)
        {}


    //- Destructor. Possibly invokes an atomic rename
    //- (preference defined during construction)
    ~OFstream();


    // Member Functions

        //- Read/write access to the name of the stream
        using OSstream::name;


    // STL stream

        //- Const access to underlying std::ostream
        virtual const std::ostream& stdStream() const override;

        //- Access to underlying std::ostream
        virtual std::ostream& stdStream() override;

        //- Rewind the stream so that it may be written again.
        //- Reopens the file (truncation)
        virtual void rewind();


    // Output stream modes

        //- True if opened in append mode \em and file already existed
        bool is_appending() const noexcept
        {
            return ofstreamPointer::is_appending();
        }

        //- True if file creation behaves as atomic
        bool is_atomic() const noexcept
        {
            return ofstreamPointer::is_atomic();
        }


    // Print

        //- Print stream description
        void print(Ostream& os) const override;


    // Additional constructors and methods
    #ifdef Foam_IOstream_extras

        //- Construct from pathname, format (version, compression)
        FOAM_DEPRECATED_FOR(2022-09, "Construct without specifying version")
        OFstream
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED,
            IOstreamOption::appendType append = IOstreamOption::NO_APPEND
        )
        :
            OFstream(pathname, IOstreamOption(fmt, ver, cmp), append)
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global predefined null output stream "/dev/null"
extern OFstream Snull;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

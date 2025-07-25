/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::IFstream

Description
    Input from file stream as an ISstream, normally using \c std::ifstream
    for the actual input.

SourceFiles
    IFstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IFstream_H
#define Foam_IFstream_H

#include "ISstream.H"
#include "className.H"
#include "fstreamPointer.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class IFstream Declaration
\*---------------------------------------------------------------------------*/

class IFstream
:
    private Foam::ifstreamPointer,
    public Foam::ISstream
{
public:

    //- Declare type-name (with debug switch)
    ClassName("IFstream");


    // Constructors

        //- Construct from pathname, default or specified stream options
        explicit IFstream
        (
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct from pathname and format
        IFstream
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt
        )
        :
            IFstream(pathname, IOstreamOption(fmt))
        {}


    //- Destructor
    ~IFstream() = default;


    // Static Functions

        //- Get file contents from specified file (compressed/uncompressed).
        //- Returns an empty list if the file cannot be opened.
        static DynamicList<char> readContents(const fileName& pathname);

        //- Get file contents from IFstream (assumed to be rewound)
        static DynamicList<char> readContents(IFstream& ifs);


    // Member Functions

        //- Get character(s)
        using ISstream::get;

        //- Read/write access to the name of the stream
        using ISstream::name;

        //- Return the size of the underlying file (-1 on error).
        //- This corresponds to Foam::fileSize() but with extra handling of
        //- compressed files.
        //  The return type is \c std::streamsize instead of \c off_t.
        //  \note Use sparingly since it involves a file stat()!
        std::streamsize fileSize() const;


    // STL stream

        //- Access to underlying std::istream
        virtual std::istream& stdStream() override;

        //- Const access to underlying std::istream
        virtual const std::istream& stdStream() const override;

        //- Rewind the stream so that it may be read again.
        //  Includes special handling for compressed streams.
        virtual void rewind() override;


    // Print

        //- Print stream description
        virtual void print(Ostream& os) const override;


    // Member Operators

        //- Return a non-const reference to const IFstream
        //  Needed for read-constructors where the stream argument is temporary:
        //  e.g. thing thisThing(IFstream("fileName")());
        IFstream& operator()() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

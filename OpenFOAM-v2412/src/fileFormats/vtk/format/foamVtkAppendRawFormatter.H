/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::vtk::appendRawFormatter

Description
    Appended raw binary output.

SourceFiles
    foamVtkAppendRawFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_appendRawFormatter_H
#define Foam_vtk_appendRawFormatter_H

#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                   Class vtk::appendRawFormatter Declaration
\*---------------------------------------------------------------------------*/

class appendRawFormatter
:
    public vtk::formatter
{
    // Private Data Members

        static const char* name_;
        static const char* encoding_;
        static const outputOptions opts_;

        //- The current offset for append data.
        uint64_t offset_;


protected:

    // Protected Member Functions

        //- Write
        void write(const char* s, std::streamsize n);

        //- No copy construct
        appendRawFormatter(const appendRawFormatter&) = delete;

        //- No copy assignment
        void operator=(const appendRawFormatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream
        explicit appendRawFormatter(std::ostream& os);


    //- Destructor
    virtual ~appendRawFormatter() = default;


    // Member Functions

        //- The output is APPEND_BINARY.
        virtual const vtk::outputOptions& opts() const;

        //- Output name for XML type ("append")
        virtual const char* name() const;

        //- Output name for append encoding type ("raw")
        virtual const char* encoding() const;

        //- Increase the append data offset by numbytes and sizeof(uint64_t).
        //  The additional (uint64_t) size information is consistent with
        //  writeSize()
        //
        //  \return The previous data offset
        virtual uint64_t offset(const uint64_t numbytes);

        //- Write leading size for binary output
        //  \return True - format uses this information
        virtual bool writeSize(const uint64_t numbytes);

        virtual void write(const uint8_t val);
        virtual void write(const label val);
        virtual void write(const float val);
        virtual void write(const double val);

        //- A no-op for this format
        virtual void flush();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

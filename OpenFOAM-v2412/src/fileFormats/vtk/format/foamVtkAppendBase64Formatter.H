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
    Foam::vtk::appendBase64Formatter

Description
    Appended base-64 encoded binary output.
    Uses an output filter layer to write base-64 encoded content.

SourceFiles
    foamVtkAppendBase64Formatter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_appendBase64Formatter_H
#define Foam_vtk_appendBase64Formatter_H

#include "foamVtkBase64Layer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                 Class vtk::appendBase64Formatter Declaration
\*---------------------------------------------------------------------------*/

class appendBase64Formatter
:
    public vtk::foamVtkBase64Layer
{
    // Private Data Members

        static const char* name_;
        static const outputOptions opts_;

        //- The current offset for append data.
        uint64_t offset_;


    // Private Member Functions

        //- No copy construct
        appendBase64Formatter(const appendBase64Formatter&) = delete;

        //- No copy assignment
        void operator=(const appendBase64Formatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream
        explicit appendBase64Formatter(std::ostream& os);


    //- Destructor. Closes/flushes the underlying layer.
    virtual ~appendBase64Formatter();


    // Member Functions

        //- The output is APPEND_BASE64.
        virtual const vtk::outputOptions& opts() const;

        //- Output name for XML type ("append")
        virtual const char* name() const;

        //- Increase the append data offset by numbytes and sizeof(uint64_t).
        //  The additional (uint64_t) size information is consistent with
        //  writeSize().
        //
        //  \return The previous data offset
        virtual uint64_t offset(const uint64_t numbytes);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

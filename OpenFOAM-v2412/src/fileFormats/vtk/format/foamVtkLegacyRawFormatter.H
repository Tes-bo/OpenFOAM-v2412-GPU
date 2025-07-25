/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::vtk::legacyRawFormatter

Description
    Binary output for the VTK legacy format, always written as big-endian
    and with 32-bit integers.

    This format should never be used for OpenFOAM with 64-bit label sizes.

SourceFiles
    foamVtkLegacyRawFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_legacyRawFormatter_H
#define Foam_vtk_legacyRawFormatter_H

#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                   Class vtk::legacyRawFormatter Declaration
\*---------------------------------------------------------------------------*/

class legacyRawFormatter
:
    public vtk::formatter
{
    // Private Data Members

        static const char* legacyName_;
        static const vtk::outputOptions opts_;


protected:

    // Protected Member Functions

        //- Write
        void write(const char* s, std::streamsize n);

        //- No copy construct
        legacyRawFormatter(const legacyRawFormatter&) = delete;

        //- No copy assignment
        void operator=(const legacyRawFormatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream
        explicit legacyRawFormatter(std::ostream& os);


    //- Destructor
    virtual ~legacyRawFormatter() = default;


    // Member Functions

        //- The output is LEGACY_BINARY.
        virtual const vtk::outputOptions& opts() const;

        //- Name for the legacy binary output type ("BINARY")
        virtual const char* name() const;

        //- Name for the XML append encoding (unused)
        //  Currently identical to name(), but do not rely on this.
        virtual const char* encoding() const;


        //- Write leading size - a no-op for legacy binary output
        //  \return False - never used by this format
        virtual bool writeSize(const uint64_t ignored);

        virtual void write(const uint8_t val);
        virtual void write(const label val);
        virtual void write(const float val);
        virtual void write(const double val);

        //- Write a newline to the output
        virtual void flush();


    // Disable some XML-only methods

        inline virtual formatter& endTag(vtk::fileTag) { return *this; }
        inline virtual formatter& endDataArray() { return *this; }
        inline virtual formatter& endFieldData() { return *this; }
        inline virtual formatter& endCellData() { return *this; }
        inline virtual formatter& endPointData() { return *this; }
        inline virtual formatter& endPiece() { return *this; }
        inline virtual formatter& endVTKFile() { return *this; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

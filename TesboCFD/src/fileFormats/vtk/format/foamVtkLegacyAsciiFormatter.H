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
    Foam::vtk::legacyAsciiFormatter

Description
    Formatting as per Foam::vtk::asciiFormatter, but with
    naming for legacy output.

SourceFiles
    foamVtkLegacyAsciiFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_legacyAsciiFormatter_H
#define Foam_vtk_legacyAsciiFormatter_H

#include "foamVtkAsciiFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                  Class vtk::legacyAsciiFormatter Declaration
\*---------------------------------------------------------------------------*/

class legacyAsciiFormatter
:
    public vtk::asciiFormatter
{
    // Private Data Members

        static const char* legacyName_;
        static const outputOptions opts_;


    // Private Member Functions

        //- No copy construct
        legacyAsciiFormatter(const legacyAsciiFormatter&) = delete;

        //- No copy assignment
        void operator=(const legacyAsciiFormatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream, use default precision
        explicit legacyAsciiFormatter(std::ostream& os);

        //- Construct and attach to an output stream, use specified precision
        legacyAsciiFormatter(std::ostream& os, unsigned prec);


    //- Destructor
    virtual ~legacyAsciiFormatter() = default;


    // Member Functions

        //- The output is LEGACY_ASCII.
        virtual const vtk::outputOptions& opts() const;

        //- Name for the legacy ascii output type ("ASCII")
        virtual const char* name() const;

        //- Name for the XML append encoding (unused).
        //  Currently identical to name(), but do not rely on this.
        virtual const char* encoding() const;


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

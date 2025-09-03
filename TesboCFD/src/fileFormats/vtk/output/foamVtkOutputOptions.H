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
    Foam::vtk::outputOptions

Description
    Encapsulated combinations of output format options.
    This is primarily useful when defining the output type based on some
    command-line arguments or dictionary contents.
    However, it can also be a useful alternative to using the underlying
    enumeration directly, since this class provides additional methods
    not possible with an enum.

SourceFiles
    foamVtkOutputOptions.C
    foamVtkOutputOptionsI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_outputOptions_H
#define Foam_vtk_outputOptions_H

#include "foamVtkOutput.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::outputOptions Declaration
\*---------------------------------------------------------------------------*/

class outputOptions
{
    // Private Member Data

        /// //- The output style tuning
        /// enum styleType
        /// {
        ///     NONE   = 0x00,  //!< Normal
        ///     HEADER = 0x01   //!< Emit xml header
        /// };

        //- The output format type
        vtk::formatType fmtType_;

        //- ASCII write precision
        mutable unsigned precision_;


public:

    // Constructors

        //- Default construct - XML insitu ASCII format with default precision
        inline outputOptions();

        //- Construct with specified format and default (ASCII) precision
        //
        //  \note This constructor should remain non-explicit.
        inline outputOptions(enum vtk::formatType fmtType);

        //- Construct with specified format and (ASCII) write precision
        inline outputOptions(enum vtk::formatType fmtType, unsigned prec);


    // Selectors

        //- Return new formatter based on the selected output options
        inline autoPtr<formatter> newFormatter(std::ostream& os) const;


    // Member Functions

    // Access

        //- The output format type
        inline vtk::formatType fmt() const noexcept;

        //- The file extension (legacy or xml) for the given content-type
        inline word ext(vtk::fileTag contentType) const;

        //- True if writer uses legacy file format
        inline bool legacy() const noexcept;

        //- True if writer uses XML file format (non-legacy)
        inline bool xml() const noexcept;

        //- True if output format uses an append mode
        inline bool append() const noexcept;

        //- True if output format does not use an append mode
        inline bool insitu() const noexcept;

        //- True if output format is ASCII
        inline bool ascii() const noexcept;

        //- Return the ASCII write precision
        inline unsigned precision() const noexcept;


    // Edit

        //- Toggle ASCII mode on/off.
        //  In XML append mode, this switches between base64 and raw binary.
        //  In XML inline mode, this switches between ASCII and base64.
        //  In legacy mode, this switches between ASCII and binary.
        //  \return outputOptions for chaining
        outputOptions& ascii(bool on);

        //- Toggle append mode on/off.
        //  \return outputOptions for chaining
        outputOptions& append(bool on);

        //- Toggle legacy mode on/off.
        //  \return outputOptions for chaining
        outputOptions& legacy(bool on);

        //- Set the write precision to be used for new ASCII formatters
        //  \return outputOptions for chaining
        outputOptions& precision(unsigned prec);


    // Other

        //- A text description about the output option selected
        string description() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkOutputOptionsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

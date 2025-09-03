/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    Foam::ensightFile

Description
    A variant of OFstream with specialised handling for Ensight writing
    of strings, integers and floats (ASCII and BINARY).

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightFile_H
#define Foam_ensightFile_H

#include "OFstream.H"
#include "ensightFileName.H"
#include "ensightVarName.H"
#include "DynamicList.H"
#include "IndirectListBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightFile Declaration
\*---------------------------------------------------------------------------*/

class ensightFile
:
    public OFstream
{
    // Private Data

        //- Allow undef in results
        static bool allowUndef_;

        //- Value to represent undef in results (default: 1e+37, floatVGREAT)
        static float undefValue_;

        //- Transient single-file:
        //- the original file size when opened in append mode, zero otherwise
        int64_t origFileSize_;

        //- Transient single-file:
        //- the time-step file-offsets (position after "BEGIN TIME STEP").
        //  Set on initial reading and added to by beginTimeStep().
        DynamicList<int64_t> timeStepOffsets_;


    // Private Member Functions

        //- Set the ASCII output formatting etc,
        //- and handle transient single-file timestep information
        void init();


    // Constructors

        //- Construct with file name, no ensight file naming adjustment.
        //  Created as an atomic or in append-mode (single-file format).
        //  In append-mode, attempts to parse existing time-step information.
        //  Changes APPEND_APP to APPEND_ATE (file rewriting).
        ensightFile
        (
            std::nullptr_t,                     // dispatch tag
            IOstreamOption::appendType append,  // (NO_APPEND or APPEND_ATE)
            const fileName& pathname,
            IOstreamOption::streamFormat fmt
        );


public:

    // Public Data Types

        //- Ensight uses \c float not \d double for floating-point
        typedef float floatType;


    // Static Data Members

        //- The keyword "coordinates"
        static const char* const coordinates;


    // Static Functions

        //- Return a null ensightFile
        static const ensightFile& null() noexcept
        {
            return NullObjectRef<ensightFile>();
        }


    // Generated Methods

        //- No copy construct
        ensightFile(const ensightFile&) = delete;

        //- No copy assignment
        void operator=(const ensightFile&) = delete;


    // Constructors

        //- Construct from path-name.
        //  The path-name is adjusted for valid ensight file naming.
        //  Created as an atomic or in append mode (single-file format).
        //  In append mode, attempts to parse existing time-step information.
        ensightFile
        (
            //! Append mode: NO_APPEND or APPEND_ATE
            IOstreamOption::appendType append,
            const fileName& pathname,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct from path and name.
        //  Only the name portion is adjusted for valid ensight file naming.
        //  Created as an atomic or in append mode (single-file format).
        //  In append mode, attempts to parse existing time-step information.
        ensightFile
        (
            //! Append mode: NO_APPEND or APPEND_ATE
            IOstreamOption::appendType append,
            const fileName& path,
            const fileName& name,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct from path-name.
        //  The path-name is adjusted for valid ensight file naming.
        //  Created as an atomic, non-append mode.
        explicit ensightFile
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        :
            ensightFile(IOstreamOption::NO_APPEND, pathname, fmt)
        {}

        //- Construct from path and name.
        //  Only the name portion is adjusted for valid ensight file naming.
        //  Created as an atomic, non-append mode.
        ensightFile
        (
            const fileName& path,
            const fileName& name,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        :
            ensightFile(IOstreamOption::NO_APPEND, path, name, fmt)
        {}


    //- Destructor. Commits the time-step footer information (if any)
    ~ensightFile();


    // Member Functions

    // Access

        //- Return setting for whether 'undef' values are allowed in results
        static bool allowUndef() noexcept;


    // Edit

        //- Enable/disable use of \c undef keyword and value
        static bool allowUndef(bool on) noexcept;

        //- Assign the value to represent undef in the results
        //  Returns the previous value
        //  NB: do not use values larger than floatScalarVGREAT
        static float undefValue(float value) noexcept;


    // Output

        //- Write "C Binary" string for binary files (eg, geometry/measured)
        void writeBinaryHeader();

        //- Write character/string content as "%79s" or as binary (max 80 chars)
        void writeString(const char* str, size_t len);

        //- Write C-string as "%79s" or as binary (max 80 chars)
        void writeString(const char* str);

        //- Write string as "%79s" or as binary (max 80 chars)
        void writeString(const std::string& str);

        //- Write integer value with specified width or as binary
        void writeInt(const int32_t val, const int fieldWidth);

        //- Write (narrowed) integer value with specified width or as binary
        void writeInt(const int64_t val, const int fieldWidth);

        //- Write floating-point with specified width or as binary
        void writeFloat(const float val, const int fieldWidth);

        //- Write (narrowed) floating-point with specified width or as binary
        void writeFloat(const double val, const int fieldWidth);

        //- Write undef value
        void writeUndef();


        //- Write element keyword with trailing newline,
        //- optionally with undef and the value for undefined
        virtual Ostream& writeKeyword(const keyType& key) override;

        //- Writing token does not make sense
        virtual bool write(const token&) override
        {
            NotImplemented;
            return false;
        }

        //- Writing single character does not make sense
        virtual Ostream& write(const char) override
        {
            NotImplemented;
            return *this;
        }

        //- Binary write
        virtual Ostream& write(const char* buf, std::streamsize count) override;

        //- Write C-string, uses writeString()
        virtual Ostream& write(const char* str) override;

        //- Write word, uses writeString()
        virtual Ostream& write(const word& str) override;

        //- Write string, uses writeString()
        virtual Ostream& write(const std::string& str) override;

        //- Write integer value as "%10d" or as binary
        virtual Ostream& write(const int32_t val) override;

        //- Write integer value as "%10d" or as binary (narrowed to int32_t)
        virtual Ostream& write(const int64_t val) override;

        //- Write floating-point as "%12.5e" or as binary
        virtual Ostream& write(const float val) override;

        //- Write floating-point as "%12.5e" or as binary (narrowed to float)
        virtual Ostream& write(const double val) override;

        //- Add carriage return to ascii stream
        void newline();


    // Transient single-file format

        //- Write "BEGIN TIME STEP" string and newline
        //- (for transient single-file format).
        //  \returns file position after the write
        int64_t beginTimeStep();

        //- Write "END TIME STEP" string and newline
        //- (for transient single-file format)
        //  \returns file position after the write
        int64_t endTimeStep();

        //- Transient single-file:
        //- write the time-step file-offsets as footer information.
        //  Maintains the current file position to allow manual use
        //  and seamless overwriting.
        //  \return the output file position at the start of the footer
        int64_t writeTimeStepFooter();

        //- Transient single-file:
        //- forget time-step file positions (advanced use)
        void clearTimeSteps() noexcept
        {
            timeStepOffsets_.clear();
        }

        //- Transient single-file:
        //- the current number of time steps
        label nTimes() const noexcept
        {
            return timeStepOffsets_.size();
        }

        //- Transient single-file:
        //- the current file-offsets for time steps within the file
        const UList<int64_t>& timeStepOffets() const noexcept
        {
            return timeStepOffsets_;
        }


    // Convenience Output Methods

        //- Begin a part (0-based index internally).
        void beginPart(const label index);

        //- Begin a part (0-based index internally), with a description.
        //- Only used for geometry files
        void beginPart(const label index, const std::string& description);

        //- Begin a "coordinates" block. Only used for geometry files.
        void beginCoordinates(const label nparticles);

        //- Begin a "particle coordinates" block (measured data)
        void beginParticleCoordinates(const label nparticles);

        //- Write a list of integers
        //  Adds newline after each value (ascii stream)
        void writeLabels(const UList<label>& list);

        //- Write a list of integers
        //  Adds newline after each value (ascii stream)
        template<class Addr>
        void writeLabels(const IndirectListBase<label, Addr>& list);

        //- Write a list of integers as float values
        //  Adds newline after each value (ascii stream)
        void writeList(const UList<label>& field);

        //- Write a list of floats as "%12.5e" or as binary.
        //  Adds newline after each value (ascii stream)
        void writeList(const UList<float>& field);

        //- Write a list of double as "%12.5e" or as binary.
        //- (double is narrowed to float)
        //  Adds newline after each value (ascii stream)
        void writeList(const UList<double>& field);

        //- Write an indirect list of float as "%12.5e" or as binary
        //  Adds newline after each value (ascii stream)
        template<class Addr>
        void writeList(const IndirectListBase<float, Addr>& field);

        //- Write an indirect list of double as "%12.5e" or as binary.
        //- (double is narrowed to float)
        //  Adds newline after each value (ascii stream)
        template<class Addr>
        void writeList(const IndirectListBase<double, Addr>& field);


    // Other Methods

        //- Check for any NaN in the field
        static bool hasUndef(const UList<float>& field);

        //- Check for any NaN in the field
        static bool hasUndef(const UList<double>& field);

        //- Check for any NaN in the field
        template<class Addr>
        static bool hasUndef(const IndirectListBase<float, Addr>& field);

        //- Check for any NaN in the field
        template<class Addr>
        static bool hasUndef(const IndirectListBase<double, Addr>& field);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightFileTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

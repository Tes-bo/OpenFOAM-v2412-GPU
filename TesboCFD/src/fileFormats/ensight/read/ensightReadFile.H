/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::ensightReadFile

Description
    A variant of IFstream with specialised handling for Ensight reading
    of strings, integers and floats (ASCII and BINARY).

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightReadFile_H
#define Foam_ensightReadFile_H

#include "IFstream.H"
#include "IOstream.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ensightReadFile Declaration
\*---------------------------------------------------------------------------*/

class ensightReadFile
:
    public IFstream
{
    // Private Data

        //- Transient single-file:
        //- beginning of time-step footer
        //- as read from the file
        int64_t timeStepFooterBegin_;

        //- Transient single-file:
        //- the time-step file-offsets (position after "BEGIN TIME STEP")
        //- as read from the file
        List<int64_t> timeStepOffsets_;


    // Private Member Functions

        //- Read string as "%80s" or as binary
        void readString(std::string& value);

        //- Initializes read information
        //- (optional detection of "C Binary" header)
        //- and scan for transient single-file format.
        void init(bool detectFormat);

        //- No copy construct
        ensightReadFile(const ensightReadFile&) = delete;

        //- No copy assignment
        void operator=(const ensightReadFile&) = delete;


public:

    //- Debug switch
    static int debug;


    // Constructors

        //- Construct a geometry reader, auto-detecting the "C Binary" header
        //- for binary files and skipping past it.
        explicit ensightReadFile
        (
            const fileName& pathname
        );

        //- Construct from pathname, use the specified (ascii/binary) format
        ensightReadFile
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt
        );


    //- Destructor
    ~ensightReadFile() = default;


    // Static Functions

        //- Extract time step footer information (if any).
        //  \return the begin of footer position.
        static int64_t getTimeStepFooter
        (
            IFstream& is,
            //! [out] File offsets for each time step (if any)
            List<int64_t>& offsets
        );


    // Read Functions

        //- Inherit read from Istream
        using Istream::read;

        //- Binary read
        virtual Istream& read(char* buf, std::streamsize count) override;

        //- Read string as "%80s" or as binary
        virtual Istream& read(string& value) override;

        //- Read integer as "%10d" or as binary (narrowed) int
        virtual Istream& read(label& value) override;

        //- Read floating-point as "%12.5e" or as binary
        virtual Istream& read(float& value) override;

        //- Read floating-point as "%12.5e" or as a binary (narrowed) float
        virtual Istream& read(double& value) override;

        //- Read element keyword. Currently the same as read(string)
        Istream& readKeyword(string& key);


    // Special Read Functions

        //- Component-wise reading of points/coordinates.
        //- Read all x components, y components and z components.
        void readPoints(const label nPoints, List<floatVector>& points);

        //- Component-wise reading of points/coordinates.
        //- Reads x components, y components and z components.
        void readPoints(const label nPoints, List<doubleVector>& points);

        //- Read and discard specified number of elements
        template<class Type>
        void skip(label n = 1)
        {
            Type dummy;
            for (; n > 0; --n)
            {
                this->read(dummy);
            }
        }


    // Transient single-file format

        //- Transient single-file:
        //- the position of the FILE_INDEX footer
        int64_t timeStepFooterBegin() const noexcept
        {
            return timeStepFooterBegin_;
        }

        //- Transient single-file:
        //- the number of time steps within the file
        label nTimes() const noexcept
        {
            return timeStepOffsets_.size();
        }

        //- Transient single-file:
        //- the file-offsets for time steps within the file
        const UList<int64_t>& timeStepOffets() const noexcept
        {
            return timeStepOffsets_;
        }

        //- Transient single-file:
        //- seek to the file position corresponding to the given time index.
        bool seekTime(const label timeIndex);


    // Housekeeping

        //- Detect if the file is \em binary by testing for initial
        //- "(C|Fortran) Binary"
        FOAM_DEPRECATED_FOR(2024-05, "detected on construct")
        static IOstreamOption::streamFormat
        detectBinaryHeader(const fileName& pathname)
        {
            ensightReadFile reader(pathname);
            return reader.format();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

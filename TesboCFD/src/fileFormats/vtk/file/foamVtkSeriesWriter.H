/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::vtk::seriesWriter

Description
    Provides a means of accumulating and generating VTK file series.

    The VTK file series format is a simple JSON format with the following
    type of content:
    \verbatim
    {
      "file-series-version" : "1.0",
      "files": [
        { "name" : "file1.vtk", "time" : 10 },
        { "name" : "file2.vtk", "time" : 20 },
        { "name" : "file3.vtk", "time" : 30 },
      ]
    }
    \endverbatim

    The append() operations include various sanity checks.
    Entries with an empty name are ignored.
    If an entry with an identical name already exists, its place
    will be overwritten with the new time value.

SourceFiles
    foamVtkSeriesWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_seriesWriter_H
#define Foam_vtk_seriesWriter_H

#include <fstream>
#include "foamVtkOutputOptions.H"
#include "instant.H"
#include "fileNameInstant.H"
#include "DynamicList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                       Class vtk::seriesWriter Declaration
\*---------------------------------------------------------------------------*/

class seriesWriter
{
    // Private Member Data

        //- A list of file/time entries
        DynamicList<fileNameInstant> entries_;

        //- Hash of existing (known) file names
        HashSet<fileName> existing_;

        //- Append the specified file/time instant.
        //  Overwrites existing entry that has the same name,
        //  does not append empty names.
        bool appendCheck(fileNameInstant inst);

        //- Remove duplicate filename entries. Keeping the last one seen.
        bool removeDuplicates();


public:

    // Constructors

        //- Construct an empty series
        seriesWriter() = default;

        //- Copy construct
        seriesWriter(const seriesWriter&) = default;

        //- Move construct
        seriesWriter(seriesWriter&&) = default;

        //- Copy assignment
        seriesWriter& operator=(const seriesWriter&) = default;

        //- Move assignment
        seriesWriter& operator=(seriesWriter&&) = default;


    //- Destructor
    ~seriesWriter() = default;


    // Static Member Functions

        //- Extract the base name for a file series
        //
        //  \param outputName   The name of the data output file
        //      Eg, "somefile_0001.vtk" would extract to "somefile.vtk"
        //  \param sep    The separator used between file stem and suffix.
        static fileName base(const fileName& outputName, char sep = '_');

        //- Extract the time-varying ending of files
        //
        //  \param file   The name of the file
        //      Eg, "somefile_0001.vtk" would extract to "0001"
        //  \param sep    The separator used between file stem and suffix.
        static word suffix(const fileName& file, char sep = '_');

        //- Print file series (JSON format) for specified time instances
        //
        //  \param os     The output stream
        //  \param base   The name for the series (eg, "path/file.vtk")
        //  \param series The list of suffix/value entries
        //  \param sep    The separator used between file stem and suffix.
        static Ostream& print
        (
            Ostream& os,
            const fileName& seriesName,
            const UList<instant>& series,
            const char sep = '_'
        );

        //- Write file series (JSON format) to disk, for specified instances
        //
        //  \param base   The name for the series (eg, "path/file.vtk")
        //  \param series The list of suffix/value entries
        //  \param sep    The separator used between file stem and suffix.
        static void write
        (
            const fileName& base,
            const UList<instant>& series,
            const char sep = '_'
        );

        //- Print file series (JSON format) for specified time instances.
        //  Since the VTK file series does not currently (OCT-2018) support
        //  sub-directories, these will be stripped on output.
        //
        //  \param os     The output stream
        //  \param series The list of filename/value entries
        static Ostream& print
        (
            Ostream& os,
            const UList<fileNameInstant>& series
        );

        //- Write file series (JSON format) to disk, for specified instances
        //
        //  \param seriesName  The name for the series (eg, "path/file.vtk")
        //  \param series The list of filename/value entries
        static void write
        (
            const fileName& seriesName,
            const UList<fileNameInstant>& series
        );


    // Member Functions

        //- True if there are no data sets
        inline bool empty() const noexcept;

        //- The number of data sets
        inline label size() const noexcept;


    // Content Management

        //- Clear entries
        inline void clear();

        //- Append the specified file instant
        inline bool append(const fileNameInstant& inst);

        //- Append the specified file instant
        inline bool append(fileNameInstant&& inst);

        //- Append the specified file instant.
        inline bool append(scalar timeValue, const fileName& file);

        //- Append the specified file instant.
        inline bool append(scalar timeValue, fileName&& file);

        //- Clear contents and reload by parsing the specified file.
        //
        //  \param seriesName the base name of the series to scan, without
        //      the ".series" ending.
        //  \param checkFiles verify that the files also exist
        //  \param restartTime ignore entries with a time greater/equal
        //      to the specified restart time.
        //
        //  \return the number of entries
        label load
        (
            const fileName& seriesName,
            const bool checkFiles = false,
            const scalar restartTime = ROOTVGREAT
        );

        //- Clear contents and scan directory for files.
        //
        //  The expected xml header content is a comment with the following:
        //  \verbatim
        //     <!-- ... time='3.14159' ... -->
        //  \endverbatim
        //
        //  \param seriesName the base name of the series to scan, without
        //      the ".series" ending.
        //  \param restartTime ignore entries with a time greater/equal
        //      to the specified restart time.
        //
        //  \return the number of entries
        label scan
        (
            const fileName& seriesName,
            const scalar restartTime = ROOTVGREAT
        );

        //- Remove entries that are greater_equal the time value.
        //
        //  \return True if the contents changed
        bool removeNewer(const scalar timeValue);

        //- Sort by time value and by file name
        void sort();


    // Writing

        //- Print file series as (JSON format)
        inline void print(Ostream& os) const;

        //- Write file series as (JSON format) to disk
        inline void write(const fileName& seriesName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkSeriesWriterI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

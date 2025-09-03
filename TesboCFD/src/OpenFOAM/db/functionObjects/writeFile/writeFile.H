/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::functionObjects::writeFile

Description
    Base class for writing single files from the function objects.

Usage
    \verbatim
    <dictName>
    {
        // Inherited entries
        ...

        // Optional entries
        writePrecision    <int>;
        writeToFile       <bool>;
        useUserTime       <bool>;
        updateHeader      <bool>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property        | Description                      | Type | Reqd | Deflt
      writePrecision  | Number of decimal points | int  | no  | \<system dflt\>
      writeToFile     | Produce text file output?        | bool | no  | true
      useUserTime     | Use user time (e.g. degrees)?    | bool | no  | true
      updateHeader    | Update header on mesh changes?   | bool | no  | true
    \endtable

Note
    The file header is normally updated whenver the mesh points or
    topology changes. In some cases, the function object is actually
    unaffected by these changes.
    Use the \c updateHeader flag to override the default behaviour.

SourceFiles
    writeFile.C
    writeFileTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeFile_H
#define functionObjects_writeFile_H

#include "objectRegistry.H"
#include "OFstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                 Class functionObjects::writeFile Declaration
\*---------------------------------------------------------------------------*/

class writeFile
{
protected:

    // Protected Data

        //- Reference to the region objectRegistry
        const objectRegistry& fileObr_;

        //- Prefix
        const fileName prefix_;

        //- Name of file
        word fileName_;

        //- File pointer
        autoPtr<OFstream> filePtr_;

        //- Write precision
        label writePrecision_;

        //- Flag to enable/disable writing to file
        bool writeToFile_;

        //- Flag to update the header, e.g. on mesh changes.
        //- Default is true.
        bool updateHeader_;

        //- Flag to identify whether the header has been written
        bool writtenHeader_;

        //- Flag to use the specified user time, e.g. CA deg instead
        //- of seconds.  Default = true
        bool useUserTime_;

        //- Start time value
        scalar startTime_;

        //- File extension; default = .dat
        string ext_;


    // Protected Member Functions

        //- Initialise the output stream for writing
        void initStream(Ostream& os) const;

        //- Return the base directory for output
        fileName baseFileDir() const;

        //- Return the base directory for the current time value
        fileName baseTimeDir() const;

        //- Return the full path for the supplied file name
        fileName filePath(const fileName& fName) const;

        //- Return autoPtr to a new file using file name
        //  Note: no check for if the file already exists
        virtual autoPtr<OFstream> newFile(const fileName& fName) const;

        //- Return autoPtr to a new file for a given time
        virtual autoPtr<OFstream> newFileAtTime
        (
            const word& name,
            scalar timeValue
        ) const;

        //- Return autoPtr to a new file using the simulation start time
        virtual autoPtr<OFstream> newFileAtStartTime
        (
            const word& name
        ) const;

        //- Reset internal file pointer to new file with new name
        virtual void resetFile(const word& name);

        //- Return the value width when writing to stream with optional offset
        Omanip<int> valueWidth(const label offset = 0) const;


        //- No copy assignment
        void operator=(const writeFile&) = delete;


    // Housekeeping

        //- Deprecated(2022-09) Return autoPtr to a new file for a given time
        //
        //  \deprecated(2022-09) - use newFileAtTime function
        FOAM_DEPRECATED_FOR(2022-09, "newFileAtTime function")
        virtual autoPtr<OFstream> createFile
        (
            const word& name,
            scalar timeValue
        ) const
        {
            return newFileAtTime(name, timeValue);
        }

        //- Deprecated(2022-09) Return autoPtr to a new file
        //- using the simulation start time
        //
        //  \deprecated(2022-09) - use newFileAtStartTime function
        FOAM_DEPRECATED_FOR(2022-09, "newFileAtStartTime function")
        virtual autoPtr<OFstream> createFile
        (
            const word& name
        ) const
        {
            return newFileAtStartTime(name);
        }


public:

    //- Additional characters for writing
    static label addChars;


    // Constructors

        //- Construct from objectRegistry, prefix, fileName
        writeFile
        (
            const objectRegistry& obr,
            const fileName& prefix,
            const word& name = "undefined",
            const bool writeToFile = true,
            const string& ext = ".dat"
        );

        //- Construct from objectRegistry, prefix, fileName
        //- and read options from dictionary
        writeFile
        (
            const objectRegistry& obr,
            const fileName& prefix,
            const word& name,
            const dictionary& dict,
            const bool writeToFile = true,
            const string& ext = ".dat"
        );

        //- Construct copy
        writeFile(const writeFile& wf);


    //- Destructor
    virtual ~writeFile() = default;


    // Member Functions

        //- Read
        virtual bool read(const dictionary& dict);

        //- Set extension
        virtual const string& setExt(const string& ext);

        //- Return access to the file (if only 1)
        virtual OFstream& file();

        //- Flag to allow writing to file
        virtual bool writeToFile() const;

        //- Flag to allow writing to the file
        virtual bool canWriteToFile() const;

        //- Flag to allow resetting the file
        virtual bool canResetFile() const;

        //- Flag to allow writing the header
        virtual bool canWriteHeader() const;

        //- Return width of character stream output
        virtual label charWidth() const;

        //- Write a commented string to stream
        virtual void writeCommented(Ostream& os, const string& str) const;

        //- Write a tabbed string to stream
        virtual void writeTabbed(Ostream& os, const string& str) const;

        //- Write a commented header to stream
        virtual void writeHeader(Ostream& os, const string& str) const;

        //- Write the current time to stream
        virtual void writeCurrentTime(Ostream& os) const;

        //- Write a break marker to the stream
        virtual void writeBreak(Ostream& os) const;

        //- Write a (commented) header property and value pair
        template<class Type>
        void writeHeaderValue
        (
            Ostream& os,
            const string& property,
            const Type& value
        ) const;

        //- Write a given value to stream with the space delimiter
        template<class Type>
        void writeValue
        (
            Ostream& os,
            const Type& val
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "writeFileTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

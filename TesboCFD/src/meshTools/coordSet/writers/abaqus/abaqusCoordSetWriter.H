/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::coordSetWriters::abaqusWriter

Description
    Write coordSet(s) as Abaqus point fields

    Example usage
    \verbatim
    T
    {
        type        sets;
        setFormat   abaqus;
        fields      (T);
        sets
        {
            ...
        }
    }
    \endverbatim

    Optional format options
    \verbatim
    formatOptions
    {
        abaqus
        {
            format          ascii;

            // Optional entries

            // Custom header: $ entries are substituions
            header
            (
                "** OpenFOAM abaqus output"
                "** Project $FOAM_CASE"
                "** File $FILE_NAME"
                "** $FIELD_NAME Time t=$TIME"
            );

            // Write geometry in addition to field data
            writeGeometry   yes;

            // Null value when sample value is not found
            // Default is scalar::min
            nullValue       0;

            // Insert additional time sub-directory in the output path
            // - yes : postProcessing/<fo-name>/<time>/<file>
            // - no  : postProcessing/<fo-name>/<file>
            useTimeDir      no;

            // Available when 'useTimeDir' is 'no' to disambiguate file names

            // Time base for output file names:
            // - 'time'      : <base>.inp_<field>.<time>
            // - 'iteration' : <base>.inp_<field>.<iteration>
            timeBase        iteration;

            // Optional start counters when using timeBase iteration
            writeIndex
            (
                T 1
            );


            ...
        }
    }
    \endverbatim

SourceFiles
    abaqusCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_abaqusWriter_H
#define Foam_coordSetWriters_abaqusWriter_H

#include "coordSetWriter.H"
#include "Enum.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                        Class abaqusWriter Declaration
\*---------------------------------------------------------------------------*/

class abaqusWriter
:
    public coordSetWriter
{
public:

    // Public enumerations

        //- Enumeration for time base
        enum class timeBase { time, iter };


private:

    // Private Data

        //- Time base names
        static const Enum<timeBase> timeBaseNames_;

        //- Optional user-defined header
        List<string> outputHeader_;

        //- User flag to write the geometry
        bool writeGeometry_;

        //- Null value; default = scalar::min
        scalar nullValue_;

        //- Optional override of localTimeDir
        bool useLocalTimeDir_;

        //- Optional time base when useLocalTimeDir_ = false
        timeBase timeBase_;

        //- Write index
        HashTable<label> writeIndex_;


    // Private Member Functions

        //- Helper to replace $WORD entries in str
        string replaceUserEntries
        (
            const string& str,
            const dictionary& vars
        ) const;

        //- Write the formatted keyword to the output stream
        Ostream& writeKeyword(Ostream& os, const word& keyword) const;

        //- Append time name when useLocalTimeDir_ = false
        void appendTimeName(const word& fieldName, fileName& fName) const;


        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,      //!< Name of field
            const Field<Type>& values   //!< Local field values to write
        );

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,
            const List<Field<Type>>& fieldValues
        );

        //- Write geometry to file.
        void writeGeometry(Ostream& os, label nTracks) const;


public:

    //- Runtime type information (no debug)
    TypeNameNoDebug("abaqus");


    // Constructors

        //- Default construct
        abaqusWriter();

        //- Default construct with specified options
        explicit abaqusWriter(const dictionary& options);

        //- Construct from components
        abaqusWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        abaqusWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~abaqusWriter();


    // Member Functions

        //- Characteristic output file name - information only
        virtual fileName path() const;  // override

        declareCoordSetWriterWriteMethod(label);
        declareCoordSetWriterWriteMethod(scalar);
        declareCoordSetWriterWriteMethod(vector);
        declareCoordSetWriterWriteMethod(sphericalTensor);
        declareCoordSetWriterWriteMethod(symmTensor);
        declareCoordSetWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordSetWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

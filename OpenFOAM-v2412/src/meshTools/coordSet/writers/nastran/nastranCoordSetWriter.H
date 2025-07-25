/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::coordSetWriters::nastranWriter

Description
    Write coordSet(s) as Nastran plot lines. Does not do field data.

    The formatOptions for nastran:
    \table
        Property    | Description                       | Required | Default
        format      | short / long / free               | no  | free
    \endtable

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry
    \verbatim
    rootdir
    `-- <time>
        |-- setName0.{nas}
        `-- setName1.{nas}
    \endverbatim

SourceFiles
    nastranCoordSetWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coordSetWriters_nastranWriter_H
#define Foam_coordSetWriters_nastranWriter_H

#include "coordSetWriter.H"
#include "NASCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{

/*---------------------------------------------------------------------------*\
                        Class nastranWriter Declaration
\*---------------------------------------------------------------------------*/

class nastranWriter
:
    public coordSetWriter
{
public:

    //- File field formats
    using fieldFormat = Foam::fileFormats::NASCore::fieldFormat;

    //- Output load format
    using loadFormat = Foam::fileFormats::NASCore::loadFormat;


private:

        //- Field format (width and separator)
        fieldFormat writeFormat_;

        //- Separator used for free format
        word separator_;


    // Private Member Functions

        //- Write the formatted keyword to the output stream
        Ostream& writeKeyword(Ostream& os, const word& keyword) const;

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
    TypeNameNoDebug("nastran");


    // Constructors

        //- Default construct
        nastranWriter();

        //- Default construct with specified options
        explicit nastranWriter(const dictionary& options);

        //- Construct from components
        nastranWriter
        (
            const coordSet& coords,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );

        //- Construct from components
        nastranWriter
        (
            const UPtrList<coordSet>& tracks,
            const fileName& outputPath,
            const dictionary& options = dictionary()
        );


    //- Destructor. Calls close()
    virtual ~nastranWriter();


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

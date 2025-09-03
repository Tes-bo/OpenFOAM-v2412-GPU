/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    Foam::surfaceWriters::nastranWriter

Description
    A surface writer for the Nastran file format - both surface mesh and fields

    The formatOptions for nastran:
    \table
        Property    | Description                           | Reqd | Default
        format      | Nastran format (short/long/free)      | no  | free
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
        commonGeometry | use separate geometry files        | no  | false
        PLOAD2      | Field selection (words/regex) for PLOAD2 | no  |
        PLOAD4      | Field selection (words/regex) for PLOAD4 | no  |
        fields      | Compat: Field pairs for PLOAD2/PLOAD4 | no  |
    \endtable

    For example,
    \verbatim
    formatOptions
    {
        nastran
        {
            format  free;     // format type

            scale   1000;     // [m] -> [mm]
            fieldScale
            {
               "p.*"   0.01;  // [Pa] -> [mbar]
            }

            // Specific NASTRAN load types
            PLOAD2 ( pMean );
            PLOAD4 ( "p.*" );

            // old style specification
            fields
            (
                (pMean  PLOAD2)
                (p      PLOAD4)
            );
        }
    }
    \endverbatim

    Unless otherwise specified, all fields will be treated as PLOAD2
    output. Can optionally specify PLOAD4 output using a combination
    of \c PLOAD4 (accept) and \c PLOAD2 (deny) entries. The older \c fields
    specification is also accepted and will be transcribed to
    corresponding PLOAD4, PLOAD2 entries.

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry
    \verbatim
    rootdir
    `-- <time>
        |-- surfaceName0.{nas}
        `-- surfaceName1.{nas}
    \endverbatim

    \subsection Fields
    \verbatim
    rootdir
    `-- <time>
        `-- field0
        |   |-- surfaceName0.{bdf}
        |   `-- surfaceName1.{bdf}
        `-- field1
            |-- surfaceName0.{bdf}
            `-- surfaceName1.{bdf}
    \endverbatim

Note
    Output variable scaling does not apply to integer types such as Ids.

SourceFiles
    nastranSurfaceWriter.C
    nastranSurfaceWriterImpl.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_nastranWriter_H
#define Foam_surfaceWriters_nastranWriter_H

#include "surfaceWriter.H"
#include "NASCore.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class nastranWriter Declaration
\*---------------------------------------------------------------------------*/

class nastranWriter
:
    public surfaceWriter
{
public:

    //- File field formats
    using fieldFormat = Foam::fileFormats::NASCore::fieldFormat;

    //- Output load format
    using loadFormat = Foam::fileFormats::NASCore::loadFormat;


private:

    // Private Data

        //- Field format (width and separator)
        fieldFormat writeFormat_;

        //- Use common geometry file
        bool commonGeometry_;

        //- Separator (used for free format)
        word separator_;

        //- Explicit selection for PLOAD2 output (deselects for PLOAD4)
        wordRes pload2_;

        //- Explicit selection for PLOAD4 output
        wordRes pload4_;


    // Private Member Functions

        //- Write a coordinate
        void writeCoord
        (
            Ostream& os,
            const point& p,
            const label pointId     //!< 0-based Point Id
        ) const;

        //- Write a face element (CTRIA3 or CQUAD4)
        void writeFace
        (
            Ostream& os,
            const word& faceType,
            const labelUList& facePts,
            const label elemId,     //!< 1-based Element Id
            const label propId      //!< 1-based Property Id
        ) const;

        //- Write the surface mesh geometry, tracking face decomposition
        //  Includes SHELL/MAT information
        //
        //  \param decompOffsets  begin/end offsets (size+1) into decompFaces
        //  \param decompFaces    Non tri/quad decomposed into triangles
        void writeGeometry
        (
            Ostream& os,
            const meshedSurf& surf,
            labelList& decompOffsets,
            DynamicList<face>& decompFaces
        ) const;

        //- Write the formatted keyword to the output stream
        Ostream& writeKeyword
        (
            Ostream& os,
            const word& keyword
        ) const;

        //- Write a formatted value to the output stream
        template<class Type>
        Ostream& writeValue(Ostream& os, const Type& value) const;

        //- Write a face-based value
        template<class Type>
        Ostream& writeFaceValue
        (
            Ostream& os,
            const loadFormat format,
            const Type& value,
            const label elemId      //!< 1-based Element Id
        ) const;


        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (with debug switch)
    TypeNameNoDebug("nastran");


    // Constructors

        //- Default construct. Default FREE format
        nastranWriter();

        //- Construct with some output options. Default FREE format
        explicit nastranWriter(const dictionary& options);

        //- Construct from components
        nastranWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        nastranWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~nastranWriter() = default;


    // Member Functions

        //- Format uses faceIds as part of its output
        virtual bool usesFaceIds() const // override
        {
            return true;
        }

        //- Write surface geometry to file.
        virtual fileName write(); // override

        declareSurfaceWriterWriteMethod(label);
        declareSurfaceWriterWriteMethod(scalar);
        declareSurfaceWriterWriteMethod(vector);
        declareSurfaceWriterWriteMethod(sphericalTensor);
        declareSurfaceWriterWriteMethod(symmTensor);
        declareSurfaceWriterWriteMethod(tensor);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

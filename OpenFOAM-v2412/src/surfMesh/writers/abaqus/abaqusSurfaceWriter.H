/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::surfaceWriters::abaqusWriter

Description
    A surface writer for the ABAQUS file format - both surface mesh and fields

    The formatOptions for abaqus:
    \table
        Property    | Description                           | Reqd | Default
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
        noGeometry  | Suppress geometry output (beta feature) | no  | false
    \endtable

    For example,
    \verbatim
    formatOptions
    {
        abaqus
        {
            scale   1000;     // [m] -> [mm]
            fieldLevel
            {
                p       1e5;    // Absolute -> gauge [Pa]
            }
            fieldScale
            {
                "p.*"   0.01;   // [Pa] -> [mbar]
            }
        }
    }
    \endverbatim

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry
    \verbatim
    rootdir
    `-- surfaceName0
        `-- geometry
            `-- surfaceName0_\<time\>.abq
    \endverbatim

    \subsection Fields
    \verbatim
    rootdir
    `-- surfaceName0
        `-- field0
        |   `-- surfaceName0_\<time\>.{inp}
        `-- field1
            `-- surfaceName0_\<time\>.{inp}
    \endverbatim

Note
    Output variable scaling does not apply to integer types such as Ids.

SourceFiles
    abaqusSurfaceWriter.C
    abaqusSurfaceWriterImpl.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_abaqusWriter_H
#define Foam_surfaceWriters_abaqusWriter_H

#include "surfaceWriter.H"
#include "ABAQUSCore.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                        Class abaqusWriter Declaration
\*---------------------------------------------------------------------------*/

class abaqusWriter
:
    public surfaceWriter
{
    // Private Data

        //- BETA feature
        bool noGeometry_;

        //- Output directory layouts
        enum class outputLayoutType
        {
            BY_TIME = 0,
            BY_FIELD,
        };

        //- Output directory layout
        outputLayoutType outputLayout_;


    // Private Member Functions

        //- Write face element (tri/quad)
        void writeFace
        (
            Ostream& os,
            const labelUList& f,
            const label elemId,     //!< 1-based Element Id
            const label propId,     //!< 1-based Property Id
            const bool header = true
        ) const;

        //- Write the surface mesh geometry, tracking face decomposition
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

        //- Write a face-based value
        template<class Type>
        Ostream& writeFaceValue
        (
            Ostream& os,
            const Type& value,
            const label elemId      //!< 0-based Element Id
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
    TypeNameNoDebug("abaqus");


    // Constructors

        //- Default construct
        abaqusWriter();

        //- Construct with some output options
        explicit abaqusWriter(const dictionary& options);

        //- Construct from components
        abaqusWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        abaqusWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~abaqusWriter() = default;


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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::surfaceWriters::x3dWriter

Description
    A surfaceWriter for X3D files.

    The formatOptions for x3d:
    \table
        Property    | Description                           | Required | Default
        range       | The min/max range for colour table    | no  | automatic
        colourMap   | The colour map for rendering          | no  | coolToWarm
        compression | Use file compression                  | no  | false
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
    \endtable

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName.x3d
    \endverbatim

    \subsection Fields
    \verbatim
    rootdir
    `-- timeName
        |-- <field0>_surfaceName.x3d
        `-- <field1>_surfaceName.x3d
    \endverbatim

Note
    The range is applied after any field scaling.

SourceFiles
    x3dSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_x3dWriter_H
#define Foam_surfaceWriters_x3dWriter_H

#include "surfaceWriter.H"
#include "X3DsurfaceFormatCore.H"
#include "colourTable.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                          Class x3dWriter Declaration
\*---------------------------------------------------------------------------*/

class x3dWriter
:
    public surfaceWriter,
    protected fileFormats::X3DsurfaceFormatCore
{
    // Private Data

        //- Output stream option
        IOstreamOption streamOpt_;

        //- Range of values
        //  The lower/upper limits for the colour table output
        //  Undefined means calculate from the data
        scalarMinMax range_;

        //- Selected colour table
        const colourTable* colourTablePtr_;


    // Private Member Functions

        //- Templated write operation
        template<class Type>
        fileName writeTemplate
        (
            const word& fieldName,          //!< Name of field
            const Field<Type>& localValues  //!< Local field values to write
        );


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("x3d");


    // Constructors

        //- Default construct
        x3dWriter();

        //- Construct with some output options
        explicit x3dWriter(const dictionary& options);

        //- Construct from components
        x3dWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        x3dWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~x3dWriter() = default;


    // Member Functions

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

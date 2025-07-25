/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::surfaceWriters::rawWriter

Description
    A surfaceWriter for raw output.

    The formatOptions for raw:
    \table
        Property    | Description                           | Required | Default
        compression | Use file compression                  | no  | false
        precision   | Write precision in ascii          | no  | same as IOstream
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
        fieldLevel  | Subtract field level before scaling   | no  | empty dict
        fieldScale  | Output field scaling                  | no  | empty dict
        normal      | Write face area normal in output      | no  | false
    \endtable

    For example,
    \verbatim
    formatOptions
    {
        raw
        {
            compression on;

            scale   1000;     // [m] -> [mm]
            normal  yes;
            precision 10;
            fieldScale
            {
               "p.*"   0.01;  // [Pa] -> [mbar]
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
    `-- timeName
        `-- surfaceName.{raw}
    \endverbatim

    \subsection Fields
    \verbatim
    rootdir
    `-- timeName
        |-- <field0>_surfaceName.{raw}
        `-- <field1>_surfaceName.{raw}
    \endverbatim

SourceFiles
    rawSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_rawWriter_H
#define Foam_surfaceWriters_rawWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                          Class rawWriter Declaration
\*---------------------------------------------------------------------------*/

class rawWriter
:
    public surfaceWriter
{
    // Private Data

        //- Output stream option
        IOstreamOption streamOpt_;

        //- ASCII write precision
        unsigned precision_;

        //- Output face area normal
        const bool writeNormal_;


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
    TypeNameNoDebug("raw");


    // Constructors

        //- Default construct
        rawWriter();

        //- Construct with some output options
        explicit rawWriter(const dictionary& options);

        //- Construct from components
        rawWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        rawWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~rawWriter() = default;


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

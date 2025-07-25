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
    Foam::surfaceWriters::foamWriter

Description
    A surfaceWriter for OpenFOAM surfaces

    \verbatim
    formatOptions
    {
        foam
        {
            format      ascii;
            compression true;
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

    Format options:
    \table
        Property | Description                              | Required | Default
        format   | ascii/binary                             | no  | ascii
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
        `-- surfaceName
            |-- "points"
            |-- "faceCentres"
            `-- "faces"
    \endverbatim

    \subsection Fields
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName
            |-- scalarField
            |   |-- field
            |   `-- field
            |-- vectorField
                |-- field
                `-- field
    \endverbatim

SourceFiles
    foamSurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_foamWriter_H
#define Foam_surfaceWriters_foamWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                         Class foamWriter Declaration
\*---------------------------------------------------------------------------*/

class foamWriter
:
    public surfaceWriter
{
    // Private Data

        //- Output stream option (default: ASCII, uncompressed)
        IOstreamOption streamOpt_;


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
    TypeNameNoDebug("foam");


    // Constructors

        //- Default construct
        foamWriter();

        //- Construct with some output options
        explicit foamWriter(const dictionary& options);

        //- Construct from components
        foamWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components
        foamWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~foamWriter() = default;


    // Member Functions

        //- The surface format has geometry in a separate file.
        virtual bool separateGeometry() const
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

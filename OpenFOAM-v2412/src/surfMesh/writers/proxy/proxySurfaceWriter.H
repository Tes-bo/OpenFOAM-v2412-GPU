/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::surfaceWriters::proxyWriter

Description
    A surfaceWriter that writes the geometry via the MeshedSurfaceProxy,
    but which does not support any fields.

    \section Output file locations

    The \c rootdir normally corresponds to something like
    \c postProcessing/\<name\>

    \subsection Geometry
    \verbatim
    rootdir
    `-- timeName
        `-- surfaceName.{obj|stl|..}
    \endverbatim

Note
    The formatOptions for proxy are file-type dependent.

    The following stream options are read and passed through,
    but not necessarily used by the corresponding backend.
    \table
        Property | Description                              | Required | Default
        format   | ascii/binary                             | no  | ascii
        compression | Use file compression                  | no  | false
        scale       | Output geometry scaling               | no  | 1
        transform   | Output coordinate transform           | no  |
    \endtable

SourceFiles
    proxySurfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceWriters_proxyWriter_H
#define Foam_surfaceWriters_proxyWriter_H

#include "surfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{

/*---------------------------------------------------------------------------*\
                         Class proxyWriter Declaration
\*---------------------------------------------------------------------------*/

class proxyWriter
:
    public surfaceWriter
{
    // Private Data

        //- The file extension associated with the proxy
        word fileExtension_;

        //- Output stream option
        IOstreamOption streamOpt_;

        //- Format options
        dictionary options_;


public:

    //- Declare type-name, virtual type (without debug switch)
    TypeNameNoDebug("proxy");


    // Constructors

        //- Construct for a given extension
        explicit proxyWriter(const word& fileExt);

        //- Construct for a given extension, with some output options
        proxyWriter(const word& fileExt, const dictionary& options);

        //- Construct from components, taking extension from outputPath
        proxyWriter
        (
            const meshedSurf& surf,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );

        //- Construct from components, taking extension from outputPath
        proxyWriter
        (
            const pointField& points,
            const faceList& faces,
            const fileName& outputPath,
            bool parallel = UPstream::parRun(),
            const dictionary& options = dictionary()
        );


    //- Destructor
    virtual ~proxyWriter() = default;


    // Factory Methods

        //- Optional select construct proxy writer.
        //  Return nullptr if the specified type is not supported.
        static autoPtr<surfaceWriter> TryNew(const word& writeType);

        //- Optional select construct proxy writer with extra write options.
        //  Return nullptr if the specified type is not supported.
        static autoPtr<surfaceWriter> TryNew
        (
            const word& writeType,
            const dictionary& writeOptions
        );


    // Member Functions

        //- A separate file is required for geometry.
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

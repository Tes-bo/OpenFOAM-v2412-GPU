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

\*---------------------------------------------------------------------------*/

#include "rawSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(rawWriter);
    addToRunTimeSelectionTable(surfaceWriter, rawWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, rawWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing implementation
#include "rawSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::rawWriter);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::rawWriter::rawWriter()
:
    surfaceWriter(),
    streamOpt_(),
    precision_(IOstream::defaultPrecision()),
    writeNormal_(false)
{}


Foam::surfaceWriters::rawWriter::rawWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    streamOpt_
    (
        IOstreamOption::ASCII,
        IOstreamOption::compressionEnum("compression", options)
    ),
    precision_
    (
        options.getOrDefault("precision", IOstream::defaultPrecision())
    ),
    writeNormal_(options.getOrDefault("normal", false))
{}


Foam::surfaceWriters::rawWriter::rawWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    rawWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::rawWriter::rawWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    rawWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::rawWriter::write()
{
    checkOpen();

    // Geometry:  rootdir/<TIME>/surfaceName.raw

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext("raw");

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }


    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        const pointField& points = surf.points();
        const faceList& faces = surf.faces();
        const bool withFaceNormal = (writeNormal_ && !this->isPointData());

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os(outputFile, streamOpt_);
        os.precision(precision_);

        // Header
        {
            os  << "# geometry NO_DATA " << faces.size() << nl;
            os  << "# x y z";
            if (withFaceNormal)
            {
                writeHeaderArea(os);
            }
            os  << nl;
        }

        // Write faces centres (optionally faceArea normals)
        for (const face& f : faces)
        {
            writePoint(os, f.centre(points));
            if (withFaceNormal)
            {
                os << ' ';
                writePoint(os, f.areaNormal(points));
            }
            os << nl;
        }

        os  << nl;
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

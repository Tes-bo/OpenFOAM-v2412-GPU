/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "boundaryDataSurfaceWriter.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "Time.H"
#include "pointIOField.H"
#include "primitivePatch.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(boundaryDataWriter);
    addToRunTimeSelectionTable(surfaceWriter, boundaryDataWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, boundaryDataWriter, wordDict);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter()
:
    surfaceWriter(),
    streamOpt_(),
    header_(true),
    writeNormal_(false)
{}


Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    streamOpt_
    (
        IOstreamOption::formatEnum("format", options, IOstreamOption::ASCII),
        IOstreamOption::compressionEnum("compression", options)
    ),
    header_(options.getOrDefault("header", true)),
    writeNormal_(options.getOrDefault("normal", false))
{}


Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    boundaryDataWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::boundaryDataWriter::boundaryDataWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    boundaryDataWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceWriters::boundaryDataWriter::serialWriteGeometry
(
    const regIOobject& iopts,
    const meshedSurf& surf
)
{
    const pointField& points = surf.points();
    const faceList& faces = surf.faces();

    if (verbose_)
    {
        if (this->isPointData())
        {
            Info<< "Writing points: " << iopts.objectPath() << endl;
        }
        else
        {
            Info<< "Writing face centres: " << iopts.objectPath() << endl;
        }
    }

    // Like regIOobject::writeObject without instance() adaptation
    // since this would write to e.g. 0/ instead of postProcessing/

    autoPtr<primitivePatch> ppPtr;

    {
        OFstream os(iopts.objectPath(), streamOpt_);

        if (header_)
        {
            iopts.writeHeader(os);
        }

        if (this->isPointData())
        {
            // Just like writeData, but without copying beforehand
            os << points;
        }
        else
        {
            ppPtr.reset(new primitivePatch(SubList<face>(faces), points));

            // Just like writeData, but without copying beforehand
            os << ppPtr().faceCentres();
        }

        if (header_)
        {
            IOobject::writeEndDivider(os);
        }
    }

    if (writeNormal_ && !this->isPointData())
    {
        vectorIOField iofld
        (
            IOobject
            (
                iopts.objectPath().path()/"normal",
                iopts.db(),
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );
        iofld.note() = "face data";

        OFstream os(iofld.objectPath(), streamOpt_);

        if (header_)
        {
            iofld.writeHeader(os);
        }

        os << ppPtr().faceNormals();

        if (header_)
        {
            IOobject::writeEndDivider(os);
        }
    }
}


Foam::fileName Foam::surfaceWriters::boundaryDataWriter::write()
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    // Dummy Time to use as objectRegistry
    refPtr<Time> timePtr(Time::NewGlobalTime());

    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!isDir(surfaceDir))
        {
            mkDir(surfaceDir);
        }

        // Write sample locations
        pointIOField iopts
        (
            IOobject
            (
                surfaceDir/"points",
                *timePtr,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            )
        );
        iopts.note() = (this->isPointData() ? "point data" : "face data");

        serialWriteGeometry(iopts, surf);
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::boundaryDataWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Geometry: rootdir/surfaceName/"points"
    // Field:    rootdir/surfaceName/<TIME>/field

    fileName surfaceDir = outputPath_;

    const fileName outputFile(surfaceDir/timeName()/fieldName);

    // Implicit geometry merge()
    tmp<Field<Type>> tfield = adjustField(fieldName, mergeField(localValues));

    if (verbose_)
    {
        Info<< " to " << outputFile << endl;
    }


    // Dummy Time to use as objectRegistry
    refPtr<Time> timePtr(Time::NewGlobalTime());

    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        // Write sample locations
        {
            pointIOField iopts
            (
                IOobject
                (
                    surfaceDir/"points",
                    *timePtr,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::NO_REGISTER
                )
            );
            iopts.note() = (this->isPointData() ? "point data" : "face data");

            serialWriteGeometry(iopts, surf);
        }

        // Write field
        {
            IOField<Type> iofld
            (
                IOobject
                (
                    outputFile,
                    *timePtr,
                    IOobjectOption::NO_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::NO_REGISTER
                )
            );
            iofld.note() = (this->isPointData() ? "point data" : "face data");

            OFstream os(iofld.objectPath(), streamOpt_);

            if (header_)
            {
                iofld.writeHeader(os);
            }

            // Just like writeData, but without copying beforehand
            os << tfield();

            if (header_)
            {
                IOobject::writeEndDivider(os);
            }
        }
    }

    wroteGeom_ = true;
    return surfaceDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::boundaryDataWriter);


// ************************************************************************* //

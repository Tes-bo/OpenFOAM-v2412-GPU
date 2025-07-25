/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "starcdSurfaceWriter.H"
#include "ListOps.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "MeshedSurfaceProxy.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(starcdWriter);
    addToRunTimeSelectionTable(surfaceWriter, starcdWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, starcdWriter, wordDict);
}
}

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Emit each component
    template<class Type>
    static inline void writeData(Ostream& os, const Type& val)
    {
        for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; ++cmpt)
        {
            os  << ' ' << component(val, cmpt);
        }
        os  << nl;
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::starcdWriter::starcdWriter()
:
    surfaceWriter(),
    streamOpt_()
{}


Foam::surfaceWriters::starcdWriter::starcdWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    streamOpt_
    (
        IOstreamOption::ASCII,
        IOstreamOption::compressionEnum("compression", options)
    )
{}


Foam::surfaceWriters::starcdWriter::starcdWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    starcdWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::starcdWriter::starcdWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    starcdWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::starcdWriter::write()
{
    checkOpen();

    // Geometry:  rootdir/<TIME>/surfaceName.{inp,cel,vrt}

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext("inp");

    if (verbose_)
    {
        Info<< "Writing geometry to " << outputFile << endl;
    }

    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        const labelUList& origFaceIds = surf.faceIds();

        // Face ids (if possible)
        const labelUList& elemIds =
        (
            !ListOps::found(origFaceIds, lessOp1<label>(0))
          ? origFaceIds
          : labelUList::null()
        );

        MeshedSurfaceProxy<face>
        (
            surf.points(),
            surf.faces(),
            UList<surfZone>::null(), // one zone
            labelUList::null(),      // no faceMap
            elemIds                  // face ids
        ).write
        (
            outputFile,
            "starcd",  // Canonical selection name
            streamOpt_
        );
    }

    wroteGeom_ = true;
    return outputFile;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::starcdWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    // Separate geometry
    if (!wroteGeom_)
    {
        write();
    }

    checkOpen();

    // Field:  rootdir/<TIME>/<field>_surfaceName.usr

    fileName outputFile = outputPath_.path();
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile /= timeName();
    }

    // Append <field>_surfaceName.usr
    outputFile /= fieldName + '_' + outputPath_.name();
    outputFile.ext("usr");

    // Implicit geometry merge()
    tmp<Field<Type>> tfield = adjustField(fieldName, mergeField(localValues));

    if (verbose_)
    {
        Info<< " to " << outputFile << endl;
    }


    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        const auto& values = tfield();

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os(outputFile, streamOpt_);

        const labelUList& elemIds = surf.faceIds();

        // Possible to use faceIds?
        const bool useOrigFaceIds =
        (
            elemIds.size() == values.size()
         && !ListOps::found(elemIds, lessOp1<label>(0))
        );

        label faceIndex = 0;

        // No header, just write values
        for (const Type& val : values)
        {
            const label elemId =
                (useOrigFaceIds ? elemIds[faceIndex] : faceIndex);

            os  << (elemId + 1);  // 1-based ids
            writeData(os, val);
            ++faceIndex;
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::starcdWriter);


// ************************************************************************* //

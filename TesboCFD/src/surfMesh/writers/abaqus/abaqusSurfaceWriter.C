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

\*---------------------------------------------------------------------------*/

#include "abaqusSurfaceWriter.H"
#include "ABAQUSCore.H"
#include "IOmanip.H"
#include "ListOps.H"
#include "OSspecific.H"
#include "surfaceWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceWriters
{
    defineTypeName(abaqusWriter);
    addToRunTimeSelectionTable(surfaceWriter, abaqusWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, abaqusWriter, wordDict);
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing implementation
#include "abaqusSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::abaqusWriter);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Write connectivity as CSV list
inline static void writeConnectivity
(
    Ostream& os,
    const label elemId,
    const labelUList& elem
)
{
    os  << "  " << elemId;

    for (const label vert : elem)
    {
        os << ", " << (vert + 1);
    }

    os << nl;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::surfaceWriters::abaqusWriter::writeFace
(
    Ostream& os,
    const labelUList& f,
    const label elemId,
    const label propId,
    bool header
) const
{
    // Only called with 3 or 4 points!

    if (header)
    {
        os  << "*ELEMENT, TYPE=S" << f.size();

        if (propId >= 0)
        {
            os  << ", ELSET=_" << propId;
        }

        os  << nl;
    }

    writeConnectivity(os, elemId, f);
}


void Foam::surfaceWriters::abaqusWriter::writeGeometry
(
    Ostream& os,
    const meshedSurf& surf,
    labelList& decompOffsets,
    DynamicList<face>& decompFaces
) const
{
    const pointField& points = surf.points();
    const faceList&    faces = surf.faces();
    const labelList&   zones = surf.zoneIds();
    const labelList& elemIds = surf.faceIds();

    // Possible to use faceIds?
    bool useOrigFaceIds =
    (
        elemIds.size() == faces.size()
     && !ListOps::found(elemIds, lessOp1<label>(0))
    );

    if (useOrigFaceIds)
    {
        // Not possible with on-the-fly face decomposition
        for (const auto& f : faces)
        {
            if (f.size() > 4)
            {
                useOrigFaceIds = false;
                break;
            }
        }
    }


    os  << "** Geometry" << nl;

    os  << nl
        << "**" << nl
        << "** Points" << nl
        << "**" << nl;

    fileFormats::ABAQUSCore::writePoints(os, points);


    // Write faces, with on-the-fly decomposition (triangulation)
    decompOffsets.resize(faces.size()+1);
    decompFaces.clear();

    decompOffsets[0] = 0; // The first offset is always zero

    os  << "**" << nl
        << "** Faces" << nl
        << "**" << nl;

    // Simple tracking for change of element type/set
    labelPair prevOutput(-1, -1);

    label elemId = 0;  // The element-id
    forAll(faces, facei)
    {
        const face& f = faces[facei];

        if (useOrigFaceIds)
        {
            // When available and not decomposed
            elemId = elemIds[facei];
        }

        // 1-offset for PID
        const label propId = 1 + (facei < zones.size() ? zones[facei] : 0);

        const label n = f.size();

        bool header =
            (prevOutput.first() != n || prevOutput.second() != propId);

        if (header)
        {
            // Update values
            prevOutput.first() = n;
            prevOutput.second() = propId;
        }

        if (n == 3 || n == 4)
        {
            writeFace(os, f, ++elemId, propId, header);
        }
        else
        {
            // Decompose into tris
            prevOutput.first() = 3;

            f.triangles(points, decompFaces);

            for
            (
                label decompi = decompOffsets[facei];
                decompi < decompFaces.size();
                ++decompi
            )
            {
                writeFace
                (
                    os,
                    decompFaces[decompi],
                    ++elemId,
                    propId,
                    header
                );

                header = false;
            }
        }

        // The end offset, which is the next begin offset
        decompOffsets[facei+1] = decompFaces.size();
    }

    os  << "**" << nl
        << "**" << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::abaqusWriter::abaqusWriter()
:
    surfaceWriter(),
    noGeometry_(false),
    outputLayout_(outputLayoutType::BY_FIELD)
{}


Foam::surfaceWriters::abaqusWriter::abaqusWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    noGeometry_(options.getOrDefault("noGeometry", false)),
    outputLayout_(outputLayoutType::BY_FIELD)
{}


Foam::surfaceWriters::abaqusWriter::abaqusWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    abaqusWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::abaqusWriter::abaqusWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    abaqusWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::abaqusWriter::write()
{
    checkOpen();

    // Geometry:
    // 1) rootdir/<TIME>/surfaceName.abq
    // 2) rootdir/geometry/surfaceName_<TIME>.abq

    fileName outputFile;

    switch (outputLayout_)
    {
        case outputLayoutType::BY_TIME:
        {
            outputFile = outputPath_;
            if (useTimeDir() && !timeName().empty())
            {
                // Splice in time-directory
                outputFile =
                    outputPath_.path() / timeName() / outputPath_.name();
            }
            break;
        }
        case outputLayoutType::BY_FIELD:
        {
            outputFile = outputPath_ / "geometry" / outputPath_.name();
            if (!timeName().empty())
            {
                // Append time information to file name
                outputFile += '_' + timeName();
            }
            break;
        }
    }
    outputFile.ext("abq");

    if (verbose_)
    {
        Info<< "Writing abaqus geometry to " << outputFile << endl;
    }


    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os(outputFile);

        labelList decompOffsets;
        DynamicList<face> decompFaces;

        writeGeometry(os, surf, decompOffsets, decompFaces);
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

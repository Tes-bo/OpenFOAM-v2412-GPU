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

\*---------------------------------------------------------------------------*/

#include "nastranSurfaceWriter.H"
#include "Pair.H"
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
    defineTypeName(nastranWriter);
    addToRunTimeSelectionTable(surfaceWriter, nastranWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, nastranWriter, wordDict);
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing implementation
#include "nastranSurfaceWriterImpl.C"

// Field writing methods
defineSurfaceWriterWriteFields(Foam::surfaceWriters::nastranWriter);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Ostream& Foam::surfaceWriters::nastranWriter::writeKeyword
(
    Ostream& os,
    const word& keyword
) const
{
    return fileFormats::NASCore::writeKeyword(os, keyword, writeFormat_);
}


void Foam::surfaceWriters::nastranWriter::writeCoord
(
    Ostream& os,
    const point& p,
    const label pointId
) const
{
    fileFormats::NASCore::writeCoord(os, p, pointId, writeFormat_);
}


void Foam::surfaceWriters::nastranWriter::writeFace
(
    Ostream& os,
    const word& faceType,
    const labelUList& facePts,
    const label elemId,
    const label propId
) const
{
    // Only valid surface elements are CTRIA3 and CQUAD4

    // Fixed short/long formats:
    // 1 CQUAD4
    // 2 EID  : element ID
    // 3 PID  : property element ID; default = EID   (blank)
    // 4 G1   : grid point index - requires starting index of 1
    // 5 G2   : grid point index
    // 6 G3   : grid point index
    // 7 G4   : grid point index
    // 8 onwards - not used

    // For CTRIA3 elements, cols 7 onwards are not used

    writeKeyword(os, faceType)  << separator_;

    os.setf(std::ios_base::right);

    writeValue(os, elemId)      << separator_;
    writeValue(os, propId);

    switch (writeFormat_)
    {
        case fieldFormat::SHORT :
        {
            for (const label pointi : facePts)
            {
                writeValue(os, pointi + 1);
            }

            break;
        }

        case fieldFormat::LONG :
        {
            forAll(facePts, i)
            {
                writeValue(os, facePts[i] + 1);
                if (i == 1)
                {
                    os  << nl;
                    os.unsetf(std::ios_base::right);
                    writeKeyword(os, "");
                    os.setf(std::ios_base::right);
                }
            }

            break;
        }

        case fieldFormat::FREE :
        {
            for (const label pointi : facePts)
            {
                os  << separator_;
                writeValue(os, pointi + 1);
            }

            break;
        }
    }

    os  << nl;
    os.unsetf(std::ios_base::right);
}


void Foam::surfaceWriters::nastranWriter::writeGeometry
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

    // Not possible with on-the-fly face decomposition
    if (useOrigFaceIds)
    {
        for (const auto& f : faces)
        {
            if (f.size() > 4)
            {
                useOrigFaceIds = false;
                break;
            }
        }
    }


    // Write points

    os  << '$' << nl
        << "$ Points" << nl
        << '$' << nl;

    forAll(points, pointi)
    {
        writeCoord(os, points[pointi], pointi);
    }

    // Write faces, with on-the-fly decomposition (triangulation)
    decompOffsets.resize(faces.size()+1);
    decompFaces.clear();

    decompOffsets[0] = 0; // The first offset is always zero

    os  << '$' << nl
        << "$ Faces" << nl
        << '$' << nl;

    label elemId = 0;  // The element-id
    forAll(faces, facei)
    {
        const face& f = faces[facei];

        if (useOrigFaceIds)
        {
            elemId = elemIds[facei];
        }

        // 1-offset for PID
        const label propId = 1 + (facei < zones.size() ? zones[facei] : 0);

        if (f.size() == 3)
        {
            writeFace(os, "CTRIA3", f, ++elemId, propId);
        }
        else if (f.size() == 4)
        {
            writeFace(os, "CQUAD4", f, ++elemId, propId);
        }
        else
        {
            // Decompose into tris
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
                    "CTRIA3",
                    decompFaces[decompi],
                    ++elemId,
                    propId
                );
            }
        }

        // The end offset, which is the next begin offset
        decompOffsets[facei+1] = decompFaces.size();
    }


    //
    // SHELL/MAT information
    //

    // Zone id have been used for the PID. Find unique values.

    labelList pidsUsed = labelHashSet(surf.zoneIds()).sortedToc();
    if (pidsUsed.empty())
    {
        pidsUsed.resize(1, Zero); // fallback
    }

    for (auto pid : pidsUsed)
    {
        writeKeyword(os, "PSHELL")  << separator_;
        writeValue(os, pid+1);  // 1-offset for PID

        for (label i = 0; i < 7; ++i)
        {
            // Dummy values
            os  << separator_;
            writeValue(os, 1);
        }
        os  << nl;
    }


    // Use single material ID

    const label MID = 1;

    writeKeyword(os, "MAT1")  << separator_;
    writeValue(os, MID);

    for (label i = 0; i < 7; ++i)
    {
        // Dummy values
        os  << separator_;
        writeValue(os, "");
    }
    os << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriters::nastranWriter::nastranWriter()
:
    surfaceWriter(),
    writeFormat_(fieldFormat::FREE),
    commonGeometry_(false),
    separator_(",")  // FREE format
{}


Foam::surfaceWriters::nastranWriter::nastranWriter
(
    const dictionary& options
)
:
    surfaceWriter(options),
    writeFormat_
    (
        fileFormats::NASCore::fieldFormatNames.getOrDefault
        (
            "format",
            options,
            fieldFormat::FREE
        )
    ),
    commonGeometry_(options.getOrDefault("commonGeometry", false))
{
    if (writeFormat_ == fieldFormat::FREE)
    {
        separator_ = ",";
    }

    // Explicit PLOAD2, PLOAD4 field specification
    options.readIfPresent("PLOAD2", pload2_);
    options.readIfPresent("PLOAD4", pload4_);

    // Compatibility:
    // Optional pairs of (field name => load format)
    // Could make conditional on (pload2_.empty() && pload4_.empty())

    List<Pair<word>> fieldPairs;
    options.readIfPresent("fields", fieldPairs);

    for (const Pair<word>& item : fieldPairs)
    {
        const loadFormat format =
            fileFormats::NASCore::loadFormatNames[item.second()];

        if (format == loadFormat::PLOAD2)
        {
            pload2_.push_back(item.first());
        }
        else if (format == loadFormat::PLOAD4)
        {
            pload4_.push_back(item.first());
        }
    }
}


Foam::surfaceWriters::nastranWriter::nastranWriter
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    nastranWriter(options)
{
    open(surf, outputPath, parallel);
}


Foam::surfaceWriters::nastranWriter::nastranWriter
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel,
    const dictionary& options
)
:
    nastranWriter(options)
{
    open(points, faces, outputPath, parallel);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::surfaceWriters::nastranWriter::write()
{
    checkOpen();

    // Geometry:  rootdir/<TIME>/surfaceName.nas

    fileName outputFile = outputPath_;
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile = outputPath_.path() / timeName() / outputPath_.name();
    }
    outputFile.ext("nas");

    if (verbose_)
    {
        Info<< "Writing nastran geometry to " << outputFile << endl;
    }


    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        if (!Foam::isDir(outputFile.path()))
        {
            Foam::mkDir(outputFile.path());
        }

        OFstream os(IOstreamOption::ATOMIC, outputFile);
        fileFormats::NASCore::setPrecision(os, writeFormat_);

        os  << "TITLE=OpenFOAM " << outputPath_.name() << " geometry" << nl
            << "BEGIN BULK" << nl;

        labelList decompOffsets;
        DynamicList<face> decompFaces;

        writeGeometry(os, surf, decompOffsets, decompFaces);

        os  << "ENDDATA" << nl;
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

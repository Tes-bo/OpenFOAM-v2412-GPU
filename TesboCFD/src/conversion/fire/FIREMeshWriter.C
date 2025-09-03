/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "FIREMeshWriter.H"
#include "Time.H"
#include "Map.H"
#include "OFstream.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::fileFormats::FIREMeshWriter::binary         = false;
bool Foam::fileFormats::FIREMeshWriter::compress       = false;
bool Foam::fileFormats::FIREMeshWriter::prefixBoundary = true;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Output newline in ascii mode, no-op in binary mode
inline void newline(Foam::OSstream& os)
{
    if (os.format() == Foam::IOstreamOption::ASCII)
    {
        os  << Foam::endl;
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::fileFormats::FIREMeshWriter::writeGeometry(OSstream& os) const
{
    const faceList&   faces  = mesh_.faces();
    const pointField& points = mesh_.points();
    const cellList&   cells  = mesh_.cells();

    // Points
    // ~~~~~~

    // Set the precision of the points data to 10
    os.precision(10);

    Info<< "points: " << points.size() << endl;
    putFireLabel(os, points.size());
    newline(os);

    forAll(points, ptI)
    {
        // scaling is normally 1 (ie, none)
        putFirePoint(os, scaleFactor_ * points[ptI]);
    }
    newline(os); // readability


    // Faces
    // ~~~~~
    // OPENFOAM - faces normals are outward-facing.
    // FIRE     - faces normals are inward-facing.

    Info<< "faces:  " << faces.size() << endl;
    putFireLabel(os, faces.size());
    newline(os);

    forAll(faces, faceI)
    {
        // flip face
        putFireLabels(os, faces[faceI].reverseFace());
    }
    newline(os); // readability


    // Cells
    // ~~~~~

    Info<< "cells:  " << cells.size() << endl;
    putFireLabel(os, cells.size());
    newline(os);

    forAll(cells, cellId)
    {
        putFireLabels(os, cells[cellId]);
    }
    newline(os); // readability

    return os.good();
}


bool Foam::fileFormats::FIREMeshWriter::writeSelections(OSstream& os) const
{
    label nZones = 0;
    label nPatches = 0;

    // remap name between patches and cell-zones conflicts!

    Map<word> patchNames;
    Map<word> zoneNames;

    wordHashSet usedPatchNames;
    wordHashSet usedZoneNames;

    // boundaries, skipping empty and processor patches
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        if (patch.size() && !isA<processorPolyPatch>(patch))
        {
            ++nPatches;

            const word oldName = patch.name();
            word newName;
            if (prefixBoundary)
            {
                newName = "BND_" + oldName;

                if (usedPatchNames.found(newName))
                {
                    newName = "BND_patch" + ::Foam::name(patchI);
                }
            }
            else
            {
                newName = oldName;

                if (usedPatchNames.found(newName))
                {
                    newName = "patch" + ::Foam::name(patchI);
                }
            }

            usedPatchNames.set(newName);
            patchNames.set(patchI, newName);
        }
    }


    // cellzones, skipping empty zones
    forAll(mesh_.cellZones(), zoneI)
    {
        const cellZone& cZone = mesh_.cellZones()[zoneI];
        if (cZone.size())
        {
            ++nZones;

            const word oldName = cZone.name();
            word newName = oldName;

            if (usedPatchNames.found(newName) || usedZoneNames.found(newName))
            {
                newName = "CEL_zone" + ::Foam::name(zoneI);
            }

            usedZoneNames.set(newName);
            zoneNames.set(zoneI, newName);
        }
    }


    //
    // actually write things
    //

    putFireLabel(os, (nZones + nPatches));
    newline(os);

    // do cell zones
    forAll(mesh_.cellZones(), zoneI)
    {
        const cellZone& cZone = mesh_.cellZones()[zoneI];

        if (cZone.size())
        {
            Info<< "cellZone " << zoneI
                << " (size: "  << cZone.size()
                << ") name: "  << zoneNames[zoneI] << nl;

            putFireString(os, zoneNames[zoneI]);
            putFireLabel(os, static_cast<int>(FIRECore::cellSelection));
            newline(os);

            putFireLabels(os, cZone);
            newline(os); // readability
        }
    }

    // do boundaries, skipping empty and processor patches
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];
        if (patch.size() && !isA<processorPolyPatch>(patch))
        {
            Info<< "patch " << patchI
                << " (start: " << patch.start() << " size: " << patch.size()
                << ") name: " << patchNames[patchI]
                << endl;

            putFireString(os, patchNames[patchI]);
            putFireLabel(os, static_cast<int>(FIRECore::faceSelection));
            newline(os);

            putFireLabels(os, patch.size(), patch.start());
            newline(os); // readability
        }

        newline(os); // readability
    }

    return os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::FIREMeshWriter::FIREMeshWriter
(
    const polyMesh& mesh,
    const scalar scaleFactor
)
:
    meshWriter(mesh, scaleFactor)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::FIREMeshWriter::write(const fileName& meshName) const
{
    IOstreamOption streamOpt;
    if (binary)
    {
        streamOpt.format(IOstreamOption::BINARY);
    }
    if (compress)
    {
        streamOpt.compression(IOstreamOption::COMPRESSED);
    }

    fileName baseName(meshName);
    if (baseName.empty())
    {
        baseName = meshWriter::defaultMeshName;

        const Time& t = mesh_.time();

        if
        (
            t.timeName() != "0"
         && t.timeName() != t.constant()
        )
        {
            baseName += "_" + t.timeName();
        }
    }
    else
    {
        const word ext(baseName.ext());

        if (FIRECore::file3dExtensions.found(ext))
        {
            FIRECore::fileExt3d fireFileType = FIRECore::file3dExtensions[ext];

            if (fireFileType == FIRECore::fileExt3d::POLY_ASCII)
            {
                streamOpt = IOstreamOption
                (
                    IOstreamOption::ASCII,
                    IOstreamOption::UNCOMPRESSED
                );
            }
            else if (fireFileType == FIRECore::fileExt3d::POLY_BINARY)
            {
                streamOpt = IOstreamOption
                (
                    IOstreamOption::BINARY,
                    IOstreamOption::UNCOMPRESSED
                );
            }
            else if (fireFileType == FIRECore::fileExt3d::POLY_ASCII_Z)
            {
                streamOpt = IOstreamOption
                (
                    IOstreamOption::ASCII,
                    IOstreamOption::COMPRESSED
                );
            }
            else if (fireFileType == FIRECore::fileExt3d::POLY_BINARY_Z)
            {
                streamOpt = IOstreamOption
                (
                    IOstreamOption::BINARY,
                    IOstreamOption::COMPRESSED
                );
            }
        }

        baseName.remove_ext();
    }


    // A slight hack. Cannot generate compressed files with the desired ending
    // So create and rename later
    const fileName filename = FIRECore::fireFileName
    (
        baseName,
        (
            streamOpt.format() == IOstreamOption::BINARY
          ? FIRECore::POLY_BINARY
          : FIRECore::POLY_ASCII
        )
    );

    autoPtr<OFstream> osPtr(new OFstream(filename, streamOpt));

    if (osPtr->good())
    {
        Info<< "Writing output to ";
        if (streamOpt.compression() == IOstreamOption::COMPRESSED)
        {
            // output .fpmaz instead of .fpma
            Info<< '"' << osPtr().name().c_str() << "z\"" << endl;
        }
        else
        {
            Info<< osPtr().name() << endl;
        }

        writeGeometry(osPtr());
        writeSelections(osPtr());

        osPtr.clear();    // implicitly close the file

        if (streamOpt.compression() == IOstreamOption::COMPRESSED)
        {
            // rename .fpma.gz -> .fpmaz
            // The '.gz' is automatically added by OFstream in compression mode
            Foam::mv(filename + ".gz", filename + "z");
        }
    }
    else
    {
        Info<<"could not open file for writing " << filename << endl;
        return false;
    }

    return true;
}


// ************************************************************************* //

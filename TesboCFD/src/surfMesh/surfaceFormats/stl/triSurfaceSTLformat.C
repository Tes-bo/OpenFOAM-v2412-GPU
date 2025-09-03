/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "STLCore.H"
#include "STLReader.H"
#include "Fstream.H"
#include "triSurface.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// A file-scope helper class to expose static member(s)
// This is a temporary measure and is expected to disappear in the future
struct triSurfaceSTLCore
:
    public Foam::fileFormats::STLCore
{
    using Foam::fileFormats::STLCore::writeBinaryHeader;
};


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readSTL(const fileName& filename, bool forceBinary)
{
    // Read in the values
    fileFormats::STLReader reader
    (
        filename,
        (
            forceBinary
          ? fileFormats::STLCore::BINARY
          : fileFormats::STLCore::UNKNOWN
        )
    );

    // Get the map for stitched surface points, with merge tolerance depending
    // on the input format
    labelList pointMap;
    const label nUniquePoints = reader.mergePointsMap(pointMap);

    const auto& readpts = reader.points();
    const labelList& zoneIds = reader.zoneIds();

    pointField& pointLst = storedPoints();
    List<face_type>& faceLst = storedFaces();

    // Sizing
    pointLst.setSize(nUniquePoints);
    faceLst.setSize(zoneIds.size());

    // Assign points
    forAll(readpts, pointi)
    {
        pointLst[pointMap[pointi]] = readpts[pointi];
    }

    // Assign triangles
    label pointi = 0;
    forAll(faceLst, facei)
    {
        auto& f = faceLst[facei];

        f[0] = pointMap[pointi++];
        f[1] = pointMap[pointi++];
        f[2] = pointMap[pointi++];
        f.region() = zoneIds[facei];
    }

    // Set patch name/index.
    if (reader.stlFormat() == fileFormats::STLCore::ASCII)
    {
        const List<word>& names = reader.names();

        patches_.setSize(names.size());
        forAll(patches_, patchi)
        {
            patches_[patchi] = geometricSurfacePatch(names[patchi], patchi);
        }
    }

    return true;
}


void Foam::triSurface::writeSTLASCII
(
    const fileName& filename,
    const bool sort
) const
{
    OFstream os(filename);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }

    labelList faceMap;
    surfacePatchList patches(calcPatches(faceMap));

    if (sort)
    {
        label faceIndex = 0;
        forAll(patches, patchi)
        {
            // Print all faces belonging to this region
            const surfacePatch& patch = patches[patchi];

            os  << "solid " << patch.name() << endl;

            for
            (
                label patchFacei = 0;
                patchFacei < patch.size();
                ++patchFacei
            )
            {
                const label facei = faceMap[faceIndex++];
                const labelledTri& f = (*this)[facei];

                // Write ASCII
                STLtriangle::write
                (
                    os,
                    faceNormals()[facei],
                    points()[f[0]],
                    points()[f[1]],
                    points()[f[2]]
                );
            }

            os  << "endsolid " << patch.name() << endl;
        }

        return;
    }

    // Get patch (=compact region) per face
    labelList patchIDs(size());
    forAll(patches, patchi)
    {
        label facei = patches[patchi].start();

        forAll(patches[patchi], i)
        {
            patchIDs[faceMap[facei++]] = patchi;
        }
    }

    label currentPatchi = -1;
    forAll(*this, facei)
    {
        if (currentPatchi != patchIDs[facei])
        {
            if (currentPatchi != -1)
            {
                // Close previous solid
                os  << "endsolid " << patches[currentPatchi].name() << nl;
            }
            currentPatchi = patchIDs[facei];
            os  << "solid " << patches[currentPatchi].name() << nl;
        }

        const labelledTri& f = (*this)[facei];

        // Write ASCII
        STLtriangle::write
        (
            os,
            faceNormals()[facei],
            points()[f[0]],
            points()[f[1]],
            points()[f[2]]
        );
    }

    if (currentPatchi != -1)
    {
        os  << "endsolid " << patches[currentPatchi].name() << nl;
    }
}


void Foam::triSurface::writeSTLBINARY(const fileName& filename) const
{
    std::ofstream os(filename, std::ios::binary);

    // Write the STL header
    triSurfaceSTLCore::writeBinaryHeader(os, this->size());

    forAll(*this, facei)
    {
        const labelledTri& f = (*this)[facei];

        // Write BINARY
        STLtriangle
        (
            faceNormals()[facei],
            points()[f[0]],
            points()[f[1]],
            points()[f[2]],
            f.region()
        ).write(os);
    }
}


// ************************************************************************* //

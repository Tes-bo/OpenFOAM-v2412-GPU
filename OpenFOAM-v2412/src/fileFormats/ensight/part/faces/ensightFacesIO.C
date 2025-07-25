/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "ensightFaces.H"
#include "ensightOutput.H"
#include "InfoProxy.H"
#include "polyMesh.H"
#include "globalIndex.H"
#include "globalMeshData.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightFaces::write
(
    ensightGeoFile& os,
    const polyMesh& mesh,
    bool parallel
) const
{
    const ensightFaces& part = *this;

    parallel = parallel && Pstream::parRun();

    // Renumber the patch points/faces into unique points
    label nPoints = 0;  // Total number of points
    labelList pointToGlobal;  // local point to unique global index
    labelList uniqueMeshPointLabels;  // unique global points


    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();


    // Use the properly sorted faceIds (ensightFaces) and do NOT use
    // the faceZone or anything else directly, otherwise the
    // point-maps will not correspond.
    // - perform face-flipping later

    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(faces, part.faceIds()),
        points
    );

    if (parallel)
    {
        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                pp.meshPoints(),
                pp.meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        nPoints = globalPointsPtr().totalSize();  // nPoints (global)
    }
    else
    {
        // Non-parallel
        // - all information already available from PrimitivePatch

        nPoints = pp.meshPoints().size();
        uniqueMeshPointLabels = pp.meshPoints();

        pointToGlobal.resize_nocopy(nPoints);
        Foam::identity(pointToGlobal);
    }

    ensightOutput::Detail::writeCoordinates
    (
        os,
        part.index(),
        part.name(),
        nPoints,  // nPoints (global)
        UIndirectList<point>(points, uniqueMeshPointLabels),
        parallel  //!< Collective write?
    );


    // Renumber the faces belonging to the faceZone,
    // from local numbering to unique global index.

    faceList patchFaces(pp.localFaces());
    ListListOps::inplaceRenumber(pointToGlobal, patchFaces);

    // Also a good place to perform face flipping
    if (part.usesFlipMap())
    {
        const boolList& flip = part.flipMap();

        forAll(patchFaces, facei)
        {
            face& f = patchFaces[facei];

            if (flip[facei])
            {
                f.flip();
            }
        }
    }

    ensightOutput::writeFaceConnectivityPresorted
    (
        os,
        part,
        patchFaces,
        parallel  //!< Collective write?
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<ensightFaces>& iproxy
)
{
    const auto& part = *iproxy;

    os << part.name().c_str();

    for (label typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        os  << ' ' << ensightFaces::elemNames[etype]
            << ':' << part.total(etype);
    }
    os  << nl;

    return os;
}


// ************************************************************************* //

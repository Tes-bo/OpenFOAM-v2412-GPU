/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "extrudePatchMesh.H"

#include "createShellMesh.H"
#include "polyTopoChange.H"
#include "wallPolyPatch.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extrudePatchMesh, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extrudePatchMesh::extrudePatchMesh
(
    const word& regionName,
    const fvMesh& mesh,
    const fvPatch& p,
    const dictionary& dict
)
:
    fvMesh
    (
        IOobject
        (
            regionName,
            mesh.facesInstance(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        ),
        Zero,
        false
    ),
    extrudedPatch_(p.patch()),
    dict_(dict)
{}


Foam::extrudePatchMesh::extrudePatchMesh
(
    const fvMesh& mesh,
    const fvPatch& p,
    const dictionary& dict,
    const word& regionName,
    const polyPatchList& regionPatches
)
:
    extrudePatchMesh(regionName, mesh, p, dict)
{
    extrudeMesh(regionPatches);
}


Foam::extrudePatchMesh::extrudePatchMesh
(
    const fvMesh& mesh,
    const fvPatch& p,
    const dictionary& dict,
    const word& regionName,
    const List<polyPatch*>& regionPatches
)
:
    extrudePatchMesh(regionName, mesh, p, dict)
{
    // Acquire ownership of the pointers
    polyPatchList plist(const_cast<List<polyPatch*>&>(regionPatches));

    extrudeMesh(plist);
}


Foam::extrudePatchMesh::extrudePatchMesh
(
    const fvMesh& mesh,
    const fvPatch& p,
    const dictionary& dict,
    const word& regionName
)
:
    extrudePatchMesh(regionName, mesh, p, dict)
{
    polyPatchList regionPatches(3);
    List<dictionary> dicts(regionPatches.size());
    List<word> patchNames(regionPatches.size());
    List<word> patchTypes(regionPatches.size());

    dicts[bottomPatchID] = dict_.subDict("bottomCoeffs");
    dicts[sidePatchID] = dict_.subDict("sideCoeffs");
    dicts[topPatchID] = dict_.subDict("topCoeffs");

    forAll(dicts, patchi)
    {
        dicts[patchi].readEntry("name", patchNames[patchi]);
        dicts[patchi].readEntry("type", patchTypes[patchi]);
    }

    forAll(regionPatches, patchi)
    {
        dictionary& patchDict = dicts[patchi];
        patchDict.set("nFaces", 0);
        patchDict.set("startFace", 0);

        regionPatches.set
        (
            patchi,
            polyPatch::New
            (
                patchNames[patchi],
                patchDict,
                patchi,
                mesh.boundaryMesh()
            )
        );
    }

    extrudeMesh(regionPatches);
}


void Foam::extrudePatchMesh::extrudeMesh(const polyPatchList& regionPatches)
{
    if (this->boundaryMesh().empty())
    {
        const bool columnCells = dict_.get<bool>("columnCells");
        scalar featAngleCos = -GREAT;
        scalar featAngle = -1;
        if (dict_.readIfPresent("featureAngle", featAngle))
        {
            featAngleCos = Foam::cos(degToRad(featAngle));
        }

        bitSet nonManifoldEdge(extrudedPatch_.nEdges());
        for (label edgeI = 0; edgeI < extrudedPatch_.nInternalEdges(); edgeI++)
        {
            if (columnCells)
            {
                nonManifoldEdge.set(edgeI);
            }
            else
            {
                const auto& fcs = extrudedPatch_.edgeFaces()[edgeI];
                if (fcs.size() > 2)
                {
                    // TBD: issue #2780 : non-manifold edges get seen as
                    // internal.
                    // This bit of code can be removed once #2780 is solved.
                    nonManifoldEdge.set(edgeI);
                }
                else if (fcs.size() == 2 && featAngleCos >= -1)
                {
                    const auto& n = extrudedPatch_.faceNormals();
                    if ((n[fcs[0]] & n[fcs[1]]) < featAngleCos)
                    {
                        nonManifoldEdge.set(edgeI);
                    }
                }
            }
        }

        autoPtr<extrudeModel> model_(extrudeModel::New(dict_));

        faceList pointGlobalRegions;
        faceList pointLocalRegions;
        labelList localToGlobalRegion;

        const primitiveFacePatch pp
        (
            extrudedPatch_, extrudedPatch_.points()
        );

        createShellMesh::calcPointRegions
        (
            this->globalData(),
            pp,
            nonManifoldEdge,
            false,

            pointGlobalRegions,
            pointLocalRegions,
            localToGlobalRegion
        );


        // Per local region an originating point
        labelList localRegionPoints(localToGlobalRegion.size());
        forAll(pointLocalRegions, facei)
        {
            const face& f = extrudedPatch_.localFaces()[facei];
            const face& pRegions = pointLocalRegions[facei];
            forAll(pRegions, fp)
            {
                localRegionPoints[pRegions[fp]] = f[fp];
            }
        }

        // Calculate region normals by reducing local region normals
        pointField localRegionNormals(localToGlobalRegion.size());
        {
            pointField localSum(localToGlobalRegion.size(), Zero);

            forAll(pointLocalRegions, facei)
            {
                const face& pRegions = pointLocalRegions[facei];
                forAll(pRegions, fp)
                {
                    label localRegionI = pRegions[fp];
                    localSum[localRegionI] +=
                        extrudedPatch_.faceNormals()[facei];
                }
            }

            Map<point> globalSum(2*localToGlobalRegion.size());

            forAll(localSum, localRegionI)
            {
                label globalRegionI = localToGlobalRegion[localRegionI];
                globalSum.insert(globalRegionI, localSum[localRegionI]);
            }

            // Reduce
            Pstream::mapCombineReduce(globalSum, plusEqOp<point>());

            forAll(localToGlobalRegion, localRegionI)
            {
                label globalRegionI = localToGlobalRegion[localRegionI];
                localRegionNormals[localRegionI] = globalSum[globalRegionI];
            }
            localRegionNormals /= mag(localRegionNormals);
        }


        // Per local region an extrusion direction
        vectorField firstDisp(localToGlobalRegion.size());
        forAll(firstDisp, regionI)
        {
            //const point& regionPt = regionCentres[regionI];
            const point& regionPt = extrudedPatch_.points()
            [
                extrudedPatch_.meshPoints()
                [
                    localRegionPoints[regionI]
                ]
            ];
            const vector& n = localRegionNormals[regionI];
            firstDisp[regionI] = model_()(regionPt, n, 1) - regionPt;
        }

        const label nNewPatches = regionPatches.size();

        // Extrude engine
        createShellMesh extruder
        (
            pp,
            pointLocalRegions,
            localRegionPoints
        );
        this->clearOut();
        this->removeFvBoundary();
        // Make sure the patches index the new mesh
        polyPatchList newRegionPatches(regionPatches.size());
        forAll(regionPatches, patchi)
        {
            newRegionPatches.set
            (
                patchi,
                regionPatches[patchi].clone(this->boundaryMesh())
            );
        }
        this->addFvPatches(newRegionPatches, true);


        // At this point we have a valid mesh with 3 patches and zero cells.
        // Determine:
        // - per face the top and bottom patch (topPatchID, bottomPatchID)
        // - per edge, per face on edge the side patch (edgePatches)
        labelListList edgePatches(extrudedPatch_.nEdges());
        forAll(edgePatches, edgeI)
        {
            const labelList& eFaces = extrudedPatch_.edgeFaces()[edgeI];

            if (eFaces.size() != 2 || nonManifoldEdge.test(edgeI))
            {
                edgePatches[edgeI].setSize(eFaces.size(), sidePatchID);
            }
        }

        polyTopoChange meshMod(nNewPatches);

        extruder.setRefinement
        (
            firstDisp,                              // first displacement
            model_().expansionRatio(),
            model_().nLayers(),                       // nLayers
            labelList(extrudedPatch_.size(), topPatchID),
            labelList(extrudedPatch_.size(), bottomPatchID),
            edgePatches,
            meshMod
        );

        autoPtr<mapPolyMesh> map = meshMod.changeMesh
        (
            *this,              // mesh to change
            false               // inflate
        );

        // Update numbering on extruder.
        extruder.updateMesh(map());

        this->setInstance(this->thisDb().time().constant());
        this->write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

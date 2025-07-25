/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "holeToFace.H"
#include "transform.H"
#include "faceSet.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "OBJstream.H"
//#include "fvMesh.H"
//#include "volFields.H"
//#include "surfaceFields.H"
#include "topoDistanceData.H"
#include "FaceCellWave.H"
#include "syncTools.H"

#include "edgeTopoDistanceData.H"
#include "PatchEdgeFaceWave.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(holeToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, holeToFace, word);
    addToRunTimeSelectionTable(topoSetSource, holeToFace, istream);
    addToRunTimeSelectionTable(topoSetFaceSource, holeToFace, word);
    addToRunTimeSelectionTable(topoSetFaceSource, holeToFace, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        holeToFace,
        word,
        hole
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetFaceSource,
        holeToFace,
        istream,
        hole
    );
}


Foam::topoSetSource::addToUsageTable Foam::holeToFace::usage_
(
    holeToFace::typeName,
    "\n    Usage: holeToFace <faceSet> ((x0 y0 z0) (x1 y1 z1))\n\n"
    "    Select faces disconnecting the individual regions"
    " (each indicated by a locations).\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::label Foam::holeToFace::globalCount
//(
//    const bitSet& isMasterFace,
//    const bitSet& set
//)
//{
//    if (set.size() != isMasterFace.size())
//    {
//        FatalErrorInFunction << "problem" << exit(FatalError);
//    }
//
//    label n = 0;
//    for (const label facei : set)
//    {
//        if (isMasterFace(facei))
//        {
//            n++;
//        }
//    }
//    return returnReduce(n, sumOp<label>());
//}


//void Foam::holeToFace::checkFaceSync
//(
//    const string& setName,
//    const bitSet& set
//) const
//{
//    if (set.size() != mesh_.nFaces())
//    {
//        FatalErrorInFunction<< "problem" << exit(FatalError);
//    }
//    bitSet orSet(set);
//    syncTools::syncFaceList(mesh_, orSet, orEqOp<unsigned int>());
//    bitSet andSet(set);
//    syncTools::syncFaceList(mesh_, andSet, andEqOp<unsigned int>());
//
//    forAll(orSet, facei)
//    {
//        if (orSet[facei] != andSet[facei] || orSet[facei] != set[facei])
//        {
//            FatalErrorInFunction<< "problem for set " << setName
//                << "face:" << facei
//                << " at:" << mesh_.faceCentres()[facei]
//                << " patch:" << mesh_.boundaryMesh().whichPatch(facei)
//                << " set:" << set[facei]
//                << " orSet:" << orSet[facei]
//                << " andSet:" << andSet[facei]
//                << exit(FatalError);
//        }
//    }
//}


//void Foam::holeToFace::checkFaceSync
//(
//    const string& fldName,
//    const labelList& fld
//) const
//{
//    if (fld.size() != mesh_.nFaces())
//    {
//        FatalErrorInFunction<< "problem" << exit(FatalError);
//    }
//    labelList maxFld(fld);
//    syncTools::syncFaceList(mesh_, maxFld, maxEqOp<label>());
//    labelList minFld(fld);
//    syncTools::syncFaceList(mesh_, minFld, minEqOp<label>());
//
//    forAll(maxFld, facei)
//    {
//        if (maxFld[facei] != minFld[facei] || maxFld[facei] != fld[facei])
//        {
//            FatalErrorInFunction<< "problem for field " << fldName
//                << "face:" << facei
//                << " at:" << mesh_.faceCentres()[facei]
//                << " patch:" << mesh_.boundaryMesh().whichPatch(facei)
//                << " fld:" << fld[facei]
//                << " maxFld:" << maxFld[facei]
//                << " minFld:" << minFld[facei]
//                << exit(FatalError);
//        }
//    }
//}


//void Foam::holeToFace::checkFaceSync
//(
//    const string& fldName,
//    const List<unsigned int>& fld
//) const
//{
//    if (fld.size() != mesh_.nFaces())
//    {
//        FatalErrorInFunction<< "problem" << exit(FatalError);
//    }
//    List<unsigned int> orFld(fld);
//    syncTools::syncFaceList(mesh_, orFld, bitOrEqOp<unsigned int>());
//    List<unsigned int> andFld(fld);
//    syncTools::syncFaceList(mesh_, andFld, bitAndEqOp<unsigned int>());
//    forAll(orFld, facei)
//    {
//        if (orFld[facei] != andFld[facei] || orFld[facei] != fld[facei])
//        {
//            FatalErrorInFunction<< "problem for field " << fldName
//                << "face:" << facei
//                << " at:" << mesh_.faceCentres()[facei]
//                << " patch:" << mesh_.boundaryMesh().whichPatch(facei)
//                << " fld:" << fld[facei]
//                << " orFld:" << orFld[facei]
//                << " andFld:" << andFld[facei]
//                << exit(FatalError);
//        }
//    }
//}


//void Foam::holeToFace::writeCellField
//(
//    const word& name,
//    const labelList& labelFld
//) const
//{
//    Pout<< "Writing field " << name << endl;
//    if (labelFld.size() != mesh_.nCells())
//    {
//        FatalErrorInFunction << exit(FatalError);
//    }
//
//    const fvMesh& fvm = dynamic_cast<const fvMesh&>(mesh_);
//
//    volScalarField fld
//    (
//        IOobject
//        (
//            name,
//            fvm.time().timeName(),
//            fvm,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            IOobject::NO_REGISTER
//        ),
//        fvm,
//        dimensionedScalar(dimless, Zero)
//    );
//    forAll(labelFld, i)
//    {
//        fld[i] = labelFld[i];
//    }
//    fld.correctBoundaryConditions();
//    fld.write();
//}


//void Foam::holeToFace::writeFaceField
//(
//    const word& name,
//    const labelList& labelFld
//) const
//{
//    Pout<< "Writing field " << name << endl;
//    if (labelFld.size() != mesh_.nFaces())
//    {
//        FatalErrorInFunction << exit(FatalError);
//    }
//
//    const fvMesh& fvm = dynamic_cast<const fvMesh&>(mesh_);
//
//    surfaceScalarField fld
//    (
//        IOobject
//        (
//            name,
//            fvm.time().timeName(),
//            fvm,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            IOobject::NO_REGISTER
//        ),
//        fvm,
//        dimensionedScalar(dimless, Zero)
//    );
//    for (label i = 0; i < mesh_.nInternalFaces(); i++)
//    {
//        fld[i] = labelFld[i];
//    }
//    surfaceScalarField::Boundary& bfld = fld.boundaryFieldRef();
//    forAll(bfld, patchi)
//    {
//        fvsPatchScalarField& pfld = bfld[patchi];
//        forAll(pfld, i)
//        {
//            pfld[i] = labelFld[pfld.patch().start()+i];
//        }
//    }
//    fld.write();
//
//    // Write as faceSet as well
//    const bitSet isMasterFace(syncTools::getInternalOrMasterFaces(mesh_));
//
//    faceSet set(mesh_, name, 100);
//    label nMasters = 0;
//    forAll(labelFld, facei)
//    {
//        if (labelFld[facei] >= 0)
//        {
//            set.insert(facei);
//            if (isMasterFace(facei))
//            {
//                nMasters++;
//            }
//        }
//    }
//    Pout<< "Writing " << returnReduce(nMasters, sumOp<label>())
//        << " >= 0 faces to faceSet " << set.name() << endl;
//    set.write();
//}


void Foam::holeToFace::writeFaces
(
    const word& name,
    const bitSet& faceFld
) const
{
    mkDir(mesh_.time().timePath());
    OBJstream str(mesh_.time().timePath()/name);
    Pout<< "Writing " << faceFld.count() << " faces to " << str.name() << endl;

    for (const label facei : faceFld)
    {
        str.write(mesh_.faces()[facei], mesh_.points(), false);
    }
}


void Foam::holeToFace::calculateDistance
(
    const labelList& seedFaces,
    const bitSet& isBlockedCell,
    const bitSet& isBlockedFace,
    labelList& cellDist,
    labelList& faceDist
) const
{
    if (isBlockedCell.size() != mesh_.nCells())
    {
        FatalErrorInFunction << "Problem" << exit(FatalError);
    }
    if (isBlockedFace.size() != mesh_.nFaces())
    {
        FatalErrorInFunction << "Problem" << exit(FatalError);
    }

    //const bitSet isMasterFace(syncTools::getInternalOrMasterFaces(mesh_));

    // Field on cells and faces.
    List<topoDistanceData<label>> cellData(mesh_.nCells());
    List<topoDistanceData<label>> faceData(mesh_.nFaces());

    // Start of changes
    List<topoDistanceData<label>> seedData
    (
        seedFaces.size(),
        topoDistanceData<label>(0, 123)
    );
    //Pout<< "Seeded "
    //    << globalCount(isMasterFace, bitSet(mesh_.nFaces(), seedFaces))
    //    << " out of " << returnReduce(mesh_.nFaces(), sumOp<label>()) << endl;

    // Make sure we don't walk through inactive cells
    //Pout<< "blocking "
    //    << returnReduce(isBlockedCell.count(), sumOp<label>())
    //    << " cells" << endl;
    for (const label celli : isBlockedCell)
    {
        cellData[celli] = topoDistanceData<label>(0, 0);
    }
    //Pout<< "blocking "
    //    << globalCount(isMasterFace, isBlockedFace) << " faces" << endl;
    for (const label facei : isBlockedFace)
    {
        faceData[facei] = topoDistanceData<label>(0, 0);
    }

    // Propagate information inwards
    FaceCellWave<topoDistanceData<label>> deltaCalc
    (
        mesh_,
        seedFaces,
        seedData,
        faceData,
        cellData,
        mesh_.globalData().nTotalCells()+1
    );

    // And extract
    //bool haveWarned = false;
    forAll(cellData, celli)
    {
        if (!isBlockedCell[celli])
        {
            if (!cellData[celli].valid(deltaCalc.data()))
            {
                //if (!haveWarned)
                //{
                //    WarningInFunction
                //        << "Did not visit some cells, e.g. cell " << celli
                //        << " at " << mesh_.cellCentres()[celli] << endl;
                //    haveWarned = true;
                //}
            }
            else
            {
                cellDist[celli] = cellData[celli].distance();
            }
        }
    }

    forAll(faceDist, facei)
    {
        if (!isBlockedFace[facei])
        {
            if (!faceData[facei].valid(deltaCalc.data()))
            {
                //if (!haveWarned)
                //{
                //    WarningInFunction
                //        << "Did not visit some faces, e.g. face " << facei
                //        << " at " << mesh_.faceCentres()[facei] << endl;
                //    haveWarned = true;
                //}
            }
            else
            {
                faceDist[facei] = faceData[facei].distance();
            }
        }
    }
}


Foam::bitSet Foam::holeToFace::frontFaces
(
    const bitSet& isSurfaceFace,
    const List<unsigned int>& locationFaces,
    const bitSet& isHoleCell
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    bitSet isFrontFace(mesh_.nFaces());
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (!isSurfaceFace[facei])
        {
            const label ownHole = isHoleCell[faceOwner[facei]];
            const label neiHole = isHoleCell[faceNeighbour[facei]];

            if (ownHole != neiHole)
            {
                unsigned int masks = locationFaces[facei];
                if (masks == 0u)
                {
                    FatalErrorInFunction << "face:" << facei
                        << " at:" << mesh_.faceCentres()[facei]
                        << " not any front" << exit(FatalError);
                }

                // Count number of bits set
                const label nSet = BitOps::bit_count(masks);

                if (nSet == 1)
                {
                    isFrontFace.set(facei);
                }
            }
        }
    }

    // Get neighbouring cell data
    bitSet isHoleNeiCell(mesh_.nBoundaryFaces());
    {
        for (const polyPatch& pp : pbm)
        {
            label bFacei = pp.start()-mesh_.nInternalFaces();
            const labelUList& faceCells = pp.faceCells();

            for (const label celli : faceCells)
            {
                isHoleNeiCell[bFacei] = isHoleCell[celli];
                ++bFacei;
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, isHoleNeiCell);
    }


    for (label facei = mesh_.nInternalFaces(); facei < mesh_.nFaces(); facei++)
    {
        if (!isSurfaceFace[facei])
        {
            const label ownHole = isHoleCell[faceOwner[facei]];
            const label neiHole = isHoleNeiCell[facei-mesh_.nInternalFaces()];

            if (ownHole != neiHole)
            {
                unsigned int masks = locationFaces[facei];
                if (masks == 0u)
                {
                    FatalErrorInFunction << "face:" << facei
                        << " at:" << mesh_.faceCentres()[facei]
                        << " not any front" << exit(FatalError);
                }

                // Count number of bits set
                const label nSet = BitOps::bit_count(masks);

                if (nSet == 1)
                {
                    isFrontFace.set(facei);
                }
            }
        }
    }
    return isFrontFace;
}


Foam::bitSet Foam::holeToFace::findClosure
(
    const bitSet& isSurfaceFace,                // intersected faces
    const bitSet& isCandidateHoleCell,          // starting blockage
    const labelListList& locationCells          // cells per zone
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    if (zonePoints_.size() < 2)
    {
        FatalErrorInFunction << "single region only : "
            << flatOutput(zonePoints_) << exit(FatalError);
    }

    if (zonePoints_.size() > 31)
    {
        FatalErrorInFunction << "only support < 32 locations in mesh."
            << " Currently : " << flatOutput(zonePoints_) << exit(FatalError);
    }


    //const bitSet isMasterFace(syncTools::getInternalOrMasterFaces(mesh_));

    bitSet isHoleCell(isCandidateHoleCell);
    for (const labelList& zoneCells : locationCells)
    {
        for (const label celli : zoneCells)
        {
            if (celli != -1)
            {
                isHoleCell.unset(celli);
            }
        }
    }

    bitSet notHoleCell(isHoleCell);
    notHoleCell.flip();

    if (debug)
    {
        Pout<< "holeToFace::findClosure :"
            << " locationCells:" << flatOutput(locationCells) << nl
            << "holeToFace::findClosure :"
            << " initial blocked faces:" << isSurfaceFace.count()
            << " candidate closure cells:" << isHoleCell.count()
            << endl;
    }

    // Distance to surface for every cell/face inside isHoleCell
    labelList surfaceCellDist(mesh_.nCells(), -1);
    labelList surfaceNeiCellDist(mesh_.nBoundaryFaces(), -1);
    labelList surfaceFaceDist(mesh_.nFaces(), -1);
    {
        calculateDistance
        (
            isSurfaceFace.toc(),    // seed faces
            notHoleCell,            // no need to walk through non-hole cells
            bitSet(mesh_.nFaces()), // no blocked faces
            surfaceCellDist,
            surfaceFaceDist
        );
        syncTools::swapBoundaryCellList
        (
            mesh_,
            surfaceCellDist,
            surfaceNeiCellDist
        );
        //writeCellField("surfaceCellDistance0", surfaceCellDist);
        //writeFaceField("surfaceFaceDistance0", surfaceFaceDist);
        //checkFaceSync("surfaceFaceDistance0", surfaceFaceDist);
    }

    if (debug)
    {
        Pout<< "holeToFace::findClosure :"
            << " calculated topological distance to initial blocked faces."
            << " max distance:" << gMax(surfaceCellDist)
            << endl;
    }


    // Find faces reachable from locationCells. If the locationCell is inside
    // the blockage it will be only the faces of the cell. If it is outside
    // it does a full walk to find the reachable faces on the outside of
    // the blockage.
    // Note: actual distance itself does not matter - only if they have been
    // visited.
    List<unsigned int> locationFaces(mesh_.nFaces(), 0u);
    forAll(locationCells, zonei)
    {
        labelList cellDist(mesh_.nCells(), -1);
        labelList faceDist(mesh_.nFaces(), -1);

        labelList seedFaces;

        const labelList& zoneCells = locationCells[zonei];
        for (const label celli : zoneCells)
        {
            if (celli != -1)
            {
                cellDist[celli] = 0;
                const cell& cFaces = mesh_.cells()[celli];
                seedFaces.append(cFaces);
                UIndirectList<label>(faceDist, cFaces) = 0;
            }
        }

        // Extra par sync. Is this needed?
        {
            bitSet isSeedFace(mesh_.nFaces(), seedFaces);
            syncTools::syncFaceList
            (
                mesh_,
                isSeedFace,
                orEqOp<unsigned int>()
            );
            seedFaces = isSeedFace.toc();
        }


        calculateDistance
        (
            seedFaces,      // seed faces
            isHoleCell,     // do not walk through blocking cells
            isSurfaceFace,  // do not walk through surface
            cellDist,
            faceDist
        );
        //writeCellField("cellDistance" + Foam::name(zonei), cellDist);
        //writeFaceField("faceDistance" + Foam::name(zonei), faceDist);
        //checkFaceSync("faceDistance", faceDist);

        // Add all reached faces
        const unsigned int mask = (1u << zonei);
        forAll(faceDist, facei)
        {
            if (faceDist[facei] >= 0)
            {
                locationFaces[facei] |= mask;
            }
        }
    }

    if (debug)
    {
        writeFaces("isSurfaceFace.obj", isSurfaceFace);
    }
    //if (debug)
    //{
    //    bitSet isMultiBitFace(mesh_.nFaces());
    //    forAll(locationFaces, facei)
    //    {
    //        const unsigned int bits = locationFaces[facei];
    //        const label nSet = BitOps::bit_count(bits);
    //        if (nSet > 1)
    //        {
    //            isMultiBitFace.set(facei);
    //        }
    //    }
    //    writeFaces("isMultiBitFace.obj", isMultiBitFace);
    //}


    //// Check that there are some faces that connect to more than one zone
    //- NOT CORRECT!!! At this point the walking has only done the distance
    //                 outside the initial set of blocked faces. We'd have to
    //                 walk through all faces before we can determine.
    //bool haveLeak = false;
    //forAll(locationFaces, facei)
    //{
    //    if (!isSurfaceFace[facei])
    //    {
    //        const unsigned int bits = locationFaces[facei];
    //        const label nSet = BitOps::bit_count(bits);
    //        if (nSet > 1)
    //        {
    //            haveLeak = true;
    //
    //            if (debug)
    //            {
    //                // Collect points
    //                DynamicList<pointField> connected;
    //                forAll(zonePoints_, zonei)
    //                {
    //                    if (bits & (1u << zonei))
    //                    {
    //                        connected.append(zonePoints_[zonei]);
    //                    }
    //                }
    //                Pout<< "holeToFace::findClosure :"
    //                    << " found initial leak at face "
    //                    << mesh_.faceCentres()[facei]
    //                    << " between zones " << flatOutput(connected)
    //                    << endl;
    //            }
    //            break;
    //        }
    //    }
    //}
    //
    //if (!returnReduceOr(haveLeak))
    //{
    //    if (debug)
    //    {
    //        Pout<< "holeToFace::findClosure :"
    //            << " did not find leak between zones "
    //            << flatOutput(zonePoints_) << endl;
    //    }
    //    return bitSet(mesh_.nFaces());
    //}


    // Front (on outside of hole cell but not connecting multiple locations
    bitSet isFrontFace(frontFaces(isSurfaceFace, locationFaces, isHoleCell));


    // Start off eroding the cells furthest away from the surface
    label surfaceDist = gMax(surfaceCellDist);

    // Work storage
    //List<unsigned int> newLocationFaces(mesh_.nFaces());
    //bitSet newHoleCell(mesh_.nCells());

    while (surfaceDist >= 0)
    {
        // Erode cells with >= surfaceDist:
        // - unmark cell as blockage (isHoleCell)
        // - mark faces of cell as visible from inside/outside

        //newLocationFaces = locationFaces;
        //newHoleCell = isHoleCell;

        label nChanged = 0;
        for (const label facei : isFrontFace)
        {
            const label own = faceOwner[facei];
            if (isHoleCell[own])
            {
                const label ownDist = surfaceCellDist[own];
                if (ownDist >= surfaceDist)
                {
                    //newHoleCell.unset(own);
                    isHoleCell.unset(own);
                    nChanged++;

                    const cell& cFaces = mesh_.cells()[own];
                    // Set corresponding bits on faces
                    const unsigned int mask = locationFaces[facei];
                    for (const label fi : cFaces)
                    {
                        //newLocationFaces[fi] |= mask;
                        locationFaces[fi] |= mask;
                    }
                }
            }
            else if (mesh_.isInternalFace(facei))
            {
                const label nei = faceNeighbour[facei];
                const label neiDist = surfaceCellDist[nei];

                if (isHoleCell[nei] && neiDist >= surfaceDist)
                {
                    //newHoleCell[nei] = false;
                    isHoleCell[nei] = false;
                    nChanged++;

                    const cell& cFaces = mesh_.cells()[nei];
                    // Set corresponding bits on faces
                    const unsigned int mask = locationFaces[facei];
                    for (const label fi : cFaces)
                    {
                        //newLocationFaces[fi] |= mask;
                        locationFaces[fi] |= mask;
                    }
                }
            }
        }

        reduce(nChanged, sumOp<label>());


        if (debug)
        {
            Pout<< "holeToFace::findClosure :"
                << " surfaceDist:" << surfaceDist
                << " front:" << isFrontFace.count()
                << " nChangedCells:" << nChanged
                << endl;
        }


        if (nChanged == 0)
        {
            // Nothing eroded at this level. Erode cells nearer to surface.
            --surfaceDist;
        }
        else
        {
            //locationFaces = newLocationFaces;
            //isHoleCell = newHoleCell;

            // Sync locationFaces
            syncTools::syncFaceList
            (
                mesh_,
                locationFaces,
                bitOrEqOp<unsigned int>()
            );

            // Calculate new front. Never include faces that are both visible
            // from outside and inside
            isFrontFace = frontFaces(isSurfaceFace, locationFaces, isHoleCell);
        }
    }


    // Debug: dump all the end fronts
    //{
    //    forAll(locationCells, zonei)
    //    {
    //        const unsigned int mask = (1u << zonei);
    //
    //        bitSet isZoneFace(mesh_.nFaces());
    //        forAll(locationFaces, facei)
    //        {
    //            if (locationFaces[facei] & mask)
    //            {
    //                isZoneFace.set(facei);
    //            }
    //        }
    //        writeFaces("isZoneFace"+Foam::name(zonei)+".obj", isZoneFace);
    //    }
    //}



    // Find faces that are connected to more than one location
    bitSet isCommonFace(mesh_.nFaces());
    forAll(locationFaces, facei)
    {
        unsigned int masks = locationFaces[facei];
        if (masks != 0u)
        {
            // Count number of bits set
            const label nSet = BitOps::bit_count(masks);
            if (nSet >= 2)
            {
                isCommonFace.set(facei);
            }
        }
    }

    // Remove faces that are on the surface
    for (const label facei : isCommonFace)
    {
        if (surfaceFaceDist[facei] == 0)
        {
            isCommonFace.unset(facei);
        }
    }

    if (debug)
    {
        Pout<< "holeToFace::findClosure :"
            << " closure faces:" << isCommonFace.count() << endl;
    }

    return isCommonFace;
}


Foam::bitSet Foam::holeToFace::erodeSet
(
    const bitSet& isSurfaceFace,
    const bitSet& isSetFace
) const
{
    // Detect cells with lots of faces in the set. WIP. Not parallel consistent.

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    bitSet isSetCell(mesh_.nCells());
    for (const label facei : isSetFace)
    {
        isSetCell.set(faceOwner[facei]);
        if (mesh_.isInternalFace(facei))
        {
            isSetCell.set(faceNeighbour[facei]);
        }
    }

    // Count number of faces per cell. Decide if surface would improve by
    // moving set
    bitSet erodedSet(isSetFace);
    for (const label celli : isSetCell)
    {
        const cell& cFaces = mesh_.cells()[celli];

        label nBlockedFaces = 0;
        label nSurfaceFaces = 0;
        for (const label facei : cFaces)
        {
            if (erodedSet[facei])
            {
                nBlockedFaces++;
            }
            else if (isSurfaceFace[facei])
            {
                nSurfaceFaces++;
            }
        }

        if ((nSurfaceFaces + nBlockedFaces) == cFaces.size())
        {
            // Single cell already disconnected by surface intersections
            for (const label facei : cFaces)
            {
                erodedSet.unset(facei);
            }
        }
    }
    syncTools::syncFaceList
    (
        mesh_,
        erodedSet,
        andEqOp<unsigned int>()
    );
    //checkFaceSync("erodedSet", erodedSet);

    for (const label celli : isSetCell)
    {
        const cell& cFaces = mesh_.cells()[celli];

        label nBlockedFaces = 0;
        for (const label facei : cFaces)
        {
            if (erodedSet[facei])
            {
                nBlockedFaces++;
            }
        }
        if (nBlockedFaces >= cFaces.size()-2)
        {
            for (const label facei : cFaces)
            {
                erodedSet.flip(facei);
            }
        }
    }
    syncTools::syncFaceList
    (
        mesh_,
        erodedSet,
        andEqOp<unsigned int>()
    );

    if (debug)
    {
        Pout<< "holeToFace::erodeSet :"
            << " starting set:" << isSetFace.count()
            << " eroded set:" << erodedSet.count() << endl;
    }

    //checkFaceSync("erodedSet", erodedSet);
    return erodedSet;
}


void Foam::holeToFace::combine
(
    topoSet& set,
    const bitSet& isSurfaceFace,
    const bitSet& isHoleCell,
    const bool add
) const
{
    labelListList locationCells(zonePoints_.size());
    forAll(zonePoints_, zonei)
    {
        const pointField& zoneLocations = zonePoints_[zonei];
        labelList& zoneCells = locationCells[zonei];
        zoneCells.setSize(zoneLocations.size());
        forAll(zoneLocations, i)
        {
            const label celli = mesh_.findCell(zoneLocations[i]);
            zoneCells[i] = celli;

            // Check that cell has at least one unblocked face so front can
            // 'escape'.
            if (celli != -1)
            {
                const cell& cFaces = mesh_.cells()[celli];
                bool hasUnblocked = false;
                for (const label facei : cFaces)
                {
                    if (!isSurfaceFace[facei])
                    {
                        hasUnblocked = true;
                        break;
                    }
                }

                if (!hasUnblocked)
                {
                    FatalErrorInFunction
                        << "problem : location:" << zoneLocations[i]
                        << " in zone:" << zonei
                        << " is found in cell at:" << celli
                        << mesh_.cellCentres()[celli]
                        << " which is completely surrounded by blocked faces"
                        << exit(FatalError);
                }
            }
        }
    }

    bitSet isClosingFace
    (
        findClosure
        (
            isSurfaceFace,      // intersected faces
            isHoleCell,         // starting blockage
            locationCells       // cells for zonePoints_
        )
    );

    if (erode_)
    {
        isClosingFace = erodeSet(isSurfaceFace, isClosingFace);
    }

    if (debug)
    {
        writeFaces("isClosingFace.obj", isClosingFace);
        //checkFaceSync("isClosingFace", isCommonFace);
    }

    for (const label facei : isClosingFace)
    {
        addOrDelete(set, facei, add);
    }
}


Foam::List<Foam::pointField> Foam::holeToFace::expand(const pointField& pts)
{
    List<pointField> allPts(pts.size());
    forAll(pts, i)
    {
        pointField& onePt = allPts[i];
        onePt.setSize(1, pts[i]);
    }
    return allPts;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::holeToFace::holeToFace
(
    const polyMesh& mesh,
    const List<pointField>& zonePoints,
    const wordList& blockedFaceNames,
    const wordList& blockedCellNames,
    const bool erode
)
:
    topoSetFaceSource(mesh),
    zonePoints_(zonePoints),
    blockedFaceNames_(blockedFaceNames),
    blockedCellNames_(blockedCellNames),
    erode_(erode)
{}


Foam::holeToFace::holeToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetFaceSource(mesh, dict),
    zonePoints_(dict.get<List<pointField>>("points")),
    blockedFaceNames_(),
    blockedCellNames_(),
    erode_(dict.getOrDefault<bool>("erode", false))
{
    // Look for 'sets' or 'set'
    word setName;
    if (!dict.readIfPresent("faceSets", blockedFaceNames_))
    {
        if (dict.readEntry("faceSet", setName))
        {
            blockedFaceNames_.resize(1);
            blockedFaceNames_.front() = std::move(setName);
        }
    }
    if (!dict.readIfPresent("cellSets", blockedCellNames_))
    {
        if (dict.readEntry("cellSet", setName))
        {
            blockedCellNames_.resize(1);
            blockedCellNames_.front() = std::move(setName);
        }
    }
}


Foam::holeToFace::holeToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetFaceSource(mesh),
    zonePoints_(expand(pointField(is))),
    blockedFaceNames_(one{}, word(checkIs(is))),
    blockedCellNames_(),
    erode_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::holeToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    // Set of additional blocked (internal or coupled) faces
    bitSet isBlockedFace(mesh_.nFaces());
    for (const word& setName : blockedFaceNames_)
    {
        faceSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
        isBlockedFace.setMany(loadedSet.begin(), loadedSet.end());
    }

    // Optional initial blocked cells
    bitSet isCandidateCell(mesh_.nCells());
    if (blockedFaceNames_.size())
    {
        for (const word& setName : blockedCellNames_)
        {
            cellSet loadedSet(mesh_, setName, IOobject::NO_REGISTER);
            isCandidateCell.setMany(loadedSet.begin(), loadedSet.end());
        }
    }
    else
    {
        isCandidateCell = true;
    }

    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding all faces to disconnect regions: "
                << flatOutput(zonePoints_) << " ..." << endl;
        }

        combine(set, isBlockedFace, isCandidateCell, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing all faces to disconnect regions: "
                << flatOutput(zonePoints_) << " ..." << endl;
        }

        combine(set, isBlockedFace, isCandidateCell, false);
    }
}


Foam::autoPtr<Foam::mapDistribute> Foam::holeToFace::calcClosure
(
    const polyMesh& mesh,
    const List<pointField>& zonePoints,
    const labelList& blockedFaces,
    const globalIndex& globalBlockedFaces,
    const bool erode,

    labelList& closureFaces,        // local faces to close gap
    labelList& closureToBlocked
)
{
    if (blockedFaces.size() != globalBlockedFaces.localSize())
    {
        FatalErrorInFunction << "problem : blockedFaces:" << blockedFaces.size()
            << " globalBlockedFaces:" << globalBlockedFaces.localSize()
            << exit(FatalError);
    }


    // Calculate faces needed to close hole (closureFaces)
    {
        const holeToFace faceSetSource
        (
            mesh,
            zonePoints,
            wordList(0),
            wordList(0),
            erode           //false
        );
        faceSet closureFaceSet(mesh, "calcClosure", 256);

        const bitSet isBlockedFace(mesh.nFaces(), blockedFaces);
        const bitSet isActiveCell(mesh.nCells(), true);

        faceSetSource.combine
        (
            closureFaceSet,
            isBlockedFace,
            isActiveCell,
            true
        );

        closureFaces = closureFaceSet.sortedToc();
    }


    if (returnReduceAnd(closureFaces.empty()))
    {
        closureToBlocked.clear();
        return nullptr;
    }


    //- Seed edges of closureFaces patch with (global) index of blockedFace

    const indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh.faces(), closureFaces),
        mesh.points()
    );
    const edgeList& edges = pp.edges();
    const labelList& mp = pp.meshPoints();
    const label nBndEdges = pp.nEdges() - pp.nInternalEdges();

    // For all faces in blockedFaces mark the edge with a face. No special
    // handling for multiple faces sharing the edge - first one wins
    EdgeMap<label> edgeMap(pp.nEdges());
    forAll(blockedFaces, i)
    {
        const label globalBlockedi = globalBlockedFaces.toGlobal(i);
        const label facei = blockedFaces[i];
        const face& f = mesh.faces()[facei];
        forAll(f, fp)
        {
            label nextFp = f.fcIndex(fp);
            edgeMap.insert(edge(f[fp], f[nextFp]), globalBlockedi);
        }
    }
    syncTools::syncEdgeMap(mesh, edgeMap, maxEqOp<label>());



    // Seed
    DynamicList<label> initialEdges(2*nBndEdges);
    DynamicList<edgeTopoDistanceData<label, indirectPrimitivePatch>>
        initialEdgesInfo(2*nBndEdges);
    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        const edge meshE = edge(mp[e[0]], mp[e[1]]);

        auto iter = edgeMap.cfind(meshE);
        if (iter.good())
        {
            // Found edge on patch connected to blocked face. Seed with the
            // (global) index of that blocked face

            initialEdges.append(edgei);
            initialEdgesInfo.append
            (
                edgeTopoDistanceData<label, indirectPrimitivePatch>
                (
                    0,          // distance
                    iter()      // globalBlockedi
                )
            );
        }
    }

    // Data on all edges and faces
    List<edgeTopoDistanceData<label, indirectPrimitivePatch>> allEdgeInfo
    (
        pp.nEdges()
    );
    List<edgeTopoDistanceData<label, indirectPrimitivePatch>> allFaceInfo
    (
        pp.size()
    );

    // Walk
    PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        edgeTopoDistanceData<label, indirectPrimitivePatch>
    > calc
    (
        mesh,
        pp,
        initialEdges,
        initialEdgesInfo,
        allEdgeInfo,
        allFaceInfo,
        returnReduce(pp.nEdges(), sumOp<label>())+1
    );


    // Per closure face the seed face
    closureToBlocked.resize_nocopy(pp.size());
    closureToBlocked = -1;
    forAll(allFaceInfo, facei)
    {
        if (allFaceInfo[facei].valid(calc.data()))
        {
            closureToBlocked[facei] = allFaceInfo[facei].data();
        }
    }
    // Above wave only guarantees unique data on coupled edges, not on
    // coupled faces (?) so explicitly sync faces
    {
        labelList syncFld(mesh.nFaces(), -1);
        UIndirectList<label>(syncFld, pp.addressing()) = closureToBlocked;
        syncTools::syncFaceList(mesh, syncFld, maxEqOp<label>());
        closureToBlocked = UIndirectList<label>(syncFld, pp.addressing());
    }

    List<Map<label>> compactMap;
    return autoPtr<mapDistribute>::New
    (
        globalBlockedFaces,
        closureToBlocked,
        compactMap
    );
}


// ************************************************************************* //

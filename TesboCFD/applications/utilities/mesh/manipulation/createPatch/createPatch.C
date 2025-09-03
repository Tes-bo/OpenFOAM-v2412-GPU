/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

Application
    createPatch

Group
    grpMeshManipulationUtilities

Description
    Create patches out of selected boundary faces, which are either
    from existing patches or from a faceSet.

    More specifically it:
    - creates new patches (from selected boundary faces).
      Synchronise faces on coupled patches.
    - synchronises points on coupled boundaries
    - remove patches with 0 faces in them

\*---------------------------------------------------------------------------*/

#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "argList.H"
#include "Time.H"
#include "OFstream.H"
#include "meshTools.H"
#include "faceSet.H"
#include "IOPtrList.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "wordRes.H"
#include "processorMeshes.H"
#include "IOdictionary.H"
#include "regionProperties.H"
#include "faceAreaWeightAMI2D.H"
#include "fvMeshTools.H"
#include "ReadFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<dictionary>, 0);
}


word patchName(const word& name, const fvMesh& mesh0, const fvMesh& mesh1)
{
    word pName(name != "none" ? name : word::null);
    pName += mesh0.name();
    pName += "_to_";
    pName += mesh1.name();
    return pName;
}


void matchPatchFaces
(
    const word& entryName,
    const word& AMIMethod,
    const dictionary& AMIDict,
    const PtrList<fvMesh>& meshes,

    const label meshi,
    const label nSourcei,
    const label sourcei,
    const labelList& patchesi,

    const label meshj,
    const label nSourcej,
    const label sourcej,
    const labelList& patchesj,

    DynamicList<labelList>& interfaceMesh0,
    DynamicList<label>& interfaceSource0,
    DynamicList<labelList>& interfacePatch0,
    DynamicList<wordList>& interfaceNames0,
    DynamicList<List<DynamicList<label>>>& interfaceFaces0,

    DynamicList<labelList>& interfaceMesh1,
    DynamicList<label>& interfaceSource1,
    DynamicList<labelList>& interfacePatch1,
    DynamicList<wordList>& interfaceNames1,
    DynamicList<List<DynamicList<label>>>& interfaceFaces1
)
{
    // Now we have:
    // - meshi, sourcei, patchesi
    // - meshj, sourcej, patchesj


    // Attempt to match patches
    forAll(patchesi, i)
    {
        const auto& ppi = meshes[meshi].boundaryMesh()[patchesi[i]];

        forAll(patchesj, j)
        {
            const auto& ppj = meshes[meshj].boundaryMesh()[patchesj[j]];

            // Use AMI to try and find matches
            auto AMPtr(AMIInterpolation::New(AMIMethod, AMIDict));

            AMPtr->calculate(ppi, ppj, nullptr);
            if
            (
                gAverage(AMPtr->tgtWeightsSum()) > SMALL
             || gAverage(AMPtr->srcWeightsSum()) > SMALL
            )
            {
                const label inti = interfaceMesh0.size();

                Info<< "Introducing interface " << inti << " between"
                    << " mesh " << meshes[meshi].name()
                    //<< " source:" << sourcei
                    << " patch " << ppi.name()
                    << " and mesh " << meshes[meshj].name()
                    //<< " source:" << sourcej
                    << " patch " << ppj.name()
                    << endl;

                // Mesh 0
                //~~~~~~~

                auto& intMesh0 = interfaceMesh0.emplace_back(nSourcei, -1);
                intMesh0[sourcei] = meshi;

                interfaceSource0.push_back(sourcei);

                auto& intPatch0 = interfacePatch0.emplace_back(nSourcei, -1);
                intPatch0[sourcei] = ppi.index();

                auto& intNames0 = interfaceNames0.emplace_back(nSourcei);
                intNames0[sourcei] =
                    patchName(entryName, meshes[meshi], meshes[meshj]);


                // Mesh 1
                //~~~~~~~

                auto& intMesh1 = interfaceMesh1.emplace_back(nSourcej, -1);
                intMesh1[sourcej] = meshj;

                interfaceSource1.push_back(sourcej);

                auto& intPatch1 = interfacePatch1.emplace_back(nSourcej, -1);
                intPatch1[sourcej] = ppj.index();

                auto& intNames1 = interfaceNames1.emplace_back(nSourcej);
                intNames1[sourcej] =
                    patchName(entryName, meshes[meshj], meshes[meshi]);

                auto& intFaces0 = interfaceFaces0.emplace_back(nSourcei);
                DynamicList<label>& faces0 = intFaces0[sourcei];
                faces0.setCapacity(ppi.size());

                auto& intFaces1 = interfaceFaces1.emplace_back(nSourcej);
                DynamicList<label>& faces1 = intFaces1[sourcej];
                faces1.setCapacity(ppj.size());


                // Mark any mesh faces:
                // - that have a weight > 0.5.
                // - contribute to these faces (even if < 0.5)

                faces0.clear();
                faces1.clear();

                // Marked as donor of face with high weight
                scalarField targetMask;
                {
                    const scalarField& weights = AMPtr->srcWeightsSum();
                    scalarField mask(weights.size(), Zero);
                    forAll(weights, facei)
                    {
                        const scalar sum = weights[facei];
                        if (sum > 0.5)
                        {
                            mask[facei] = sum;
                        }
                    }

                    // Push source field mask to target
                    targetMask = AMPtr->interpolateToTarget(mask);
                }

                scalarField sourceMask;
                {
                    const scalarField& weights = AMPtr->tgtWeightsSum();
                    scalarField mask(weights.size(), Zero);
                    forAll(weights, facei)
                    {
                        const scalar sum = weights[facei];
                        if (sum > 0.5)
                        {
                            mask[facei] = sum;
                        }
                    }

                    // Push target mask back to source
                    sourceMask = AMPtr->interpolateToSource(mask);
                }

                {
                    const scalarField& weights = AMPtr->srcWeightsSum();
                    forAll(weights, facei)
                    {
                        if (weights[facei] > 0.5 || sourceMask[facei] > SMALL)
                        {
                            faces0.push_back(ppi.start()+facei);
                        }
                    }
                }
                {
                    const scalarField& weights = AMPtr->tgtWeightsSum();
                    forAll(weights, facei)
                    {
                        if (weights[facei] > 0.5 || targetMask[facei] > SMALL)
                        {
                            faces1.push_back(ppj.start()+facei);
                        }
                    }
                }

                faces0.shrink();
                faces1.shrink();
            }
        }
    }
}


void matchPatchFaces
(
    const bool includeOwn,
    const PtrList<fvMesh>& meshes,
    const List<DynamicList<label>>& interRegionSources,
    const List<wordList>& patchNames,
    const labelListListList& matchPatchIDs,
    const List<wordList>& AMIMethods,
    List<PtrList<dictionary>> patchInfoDicts,

    DynamicList<labelList>& interfaceMesh0,
    DynamicList<label>& interfaceSource0,
    DynamicList<labelList>& interfacePatch0,
    DynamicList<List<DynamicList<label>>>& interfaceFaces0,
    DynamicList<wordList>& interfaceNames0,

    DynamicList<labelList>& interfaceMesh1,
    DynamicList<label>& interfaceSource1,
    DynamicList<labelList>& interfacePatch1,
    DynamicList<List<DynamicList<label>>>& interfaceFaces1,
    DynamicList<wordList>& interfaceNames1
)
{
    // Add approximate matches
    forAll(meshes, meshi)
    {
        const labelListList& allPatchesi = matchPatchIDs[meshi];

        const label meshjStart(includeOwn ? meshi : meshi+1);

        for (label meshj = meshjStart; meshj < meshes.size(); meshj++)
        {
            const labelListList& allPatchesj = matchPatchIDs[meshj];

            for (const label sourcei : interRegionSources[meshi])
            {
                const labelList& patchesi = allPatchesi[sourcei];

                for (const label sourcej : interRegionSources[meshj])
                {
                    const labelList& patchesj = allPatchesj[sourcej];

                    // Now we have:
                    // - meshi, sourcei, patchesi
                    // - meshj, sourcej, patchesj

                    matchPatchFaces
                    (
                        patchNames[meshi][sourcei],
                        (
                            AMIMethods[meshi][sourcei].size()
                          ? AMIMethods[meshi][sourcei]
                          : faceAreaWeightAMI2D::typeName
                        ),
                        patchInfoDicts[meshi][sourcei],
                        meshes,

                        meshi,
                        allPatchesi.size(), // nSourcei,
                        sourcei,
                        patchesi,

                        meshj,
                        allPatchesj.size(), // nSourcej,
                        sourcej,
                        patchesj,

                        interfaceMesh0,
                        interfaceSource0,
                        interfacePatch0,
                        interfaceNames0,
                        interfaceFaces0,

                        interfaceMesh1,
                        interfaceSource1,
                        interfacePatch1,
                        interfaceNames1,
                        interfaceFaces1
                    );
                }
            }
        }
    }
}


void changePatchID
(
    const bool modify,
    const fvMesh& mesh,
    const label faceID,
    const label patchID,
    polyTopoChange& meshMod
)
{
    const label zoneID = mesh.faceZones().whichZone(faceID);

    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];

        zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
    }

    if (modify)
    {
        meshMod.setAction
        (
            polyModifyFace
            (
                mesh.faces()[faceID],               // face
                faceID,                             // face ID
                mesh.faceOwner()[faceID],           // owner
                -1,                                 // neighbour
                false,                              // flip flux
                patchID,                            // patch ID
                false,                              // remove from zone
                zoneID,                             // zone ID
                zoneFlip                            // zone flip
            )
        );
    }
    else
    {
        meshMod.setAction
        (
            polyAddFace
            (
                mesh.faces()[faceID],       // modified face
                mesh.faceOwner()[faceID],   // owner
                -1,                         // neighbour
                -1,                         // master point
                -1,                         // master edge
                faceID,                     // master face
                false,                      // face flip
                patchID,                    // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
}


void changePatchID
(
    const fvMesh& mesh,
    const labelList& faceLabels,
    const label patchID,
    bitSet& isRepatchedBoundary,
    polyTopoChange& meshMod
)
{
    for (const label facei : faceLabels)
    {
        if (mesh.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Face " << facei
                << " is not an external face of the mesh." << endl
                << "This application can only repatch"
                << " existing boundary faces." << exit(FatalError);
        }

        const bool isFirst =
            isRepatchedBoundary.set(facei-mesh.nInternalFaces());
        if (!isFirst)
        {
            static label nWarnings = 0;
            if (nWarnings == 0)
            {
                const label newPatchi = meshMod.region()[facei];
                //FatalErrorInFunction
                WarningInFunction
                    << "Face " << facei
                    << " at " << mesh.faceCentres()[facei]
                    << " marked for patch " << patchID
                    << " name " << mesh.boundaryMesh()[patchID].name()
                    << " is already marked for patch " << newPatchi
                    << " name " << mesh.boundaryMesh()[newPatchi].name()
                    << ". Creating duplicate face. Suppressing further warnings"
                    //<< exit(FatalError);
                    << endl;
            }
            nWarnings++;
        }

        changePatchID(isFirst, mesh, facei, patchID, meshMod);
    }
}


// Dump for all patches the current match
void dumpCyclicMatch(const fileName& prefix, const polyMesh& mesh)
{
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            // Dump patches
            {
                OFstream str(prefix+cycPatch.name()+".obj");
                Pout<< "Dumping " << cycPatch.name()
                    << " faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    cycPatch,
                    cycPatch.points()
                );
            }

            {
                OFstream str(prefix+nbrPatch.name()+".obj");
                Pout<< "Dumping " << nbrPatch.name()
                    << " faces to " << str.name() << endl;
                meshTools::writeOBJ
                (
                    str,
                    nbrPatch,
                    nbrPatch.points()
                );
            }


            // Lines between corresponding face centres
            OFstream str(prefix+cycPatch.name()+nbrPatch.name()+"_match.obj");
            label vertI = 0;

            Pout<< "Dumping cyclic match as lines between face centres to "
                << str.name() << endl;

            forAll(cycPatch, facei)
            {
                const point& fc0 = mesh.faceCentres()[cycPatch.start()+facei];
                meshTools::writeOBJ(str, fc0);
                vertI++;
                const point& fc1 = mesh.faceCentres()[nbrPatch.start()+facei];
                meshTools::writeOBJ(str, fc1);
                vertI++;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
}


void separateList
(
    const vectorField& separation,
    UList<vector>& field
)
{
    if (separation.size() == 1)
    {
        // Single value for all.

        forAll(field, i)
        {
            field[i] += separation[0];
        }
    }
    else if (separation.size() == field.size())
    {
        forAll(field, i)
        {
            field[i] += separation[i];
        }
    }
    else
    {
        FatalErrorInFunction
            << "Sizes of field and transformation not equal. field:"
            << field.size() << " transformation:" << separation.size()
            << abort(FatalError);
    }
}


// Synchronise points on both sides of coupled boundaries.
template<class CombineOp>
void syncPoints
(
    const polyMesh& mesh,
    pointField& points,
    const CombineOp& cop,
    const point& nullValue
)
{
    if (points.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << "Number of values " << points.size()
            << " is not equal to the number of points in the mesh "
            << mesh.nPoints() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Is there any coupled patch with transformation?
    bool hasTransformation = false;

    if (UPstream::parRun())
    {
        const labelList& procPatches = mesh.globalData().processorPatches();

        // Send
        for (const label patchi : procPatches)
        {
            const polyPatch& pp = patches[patchi];
            const auto& procPatch = refCast<const processorPolyPatch>(pp);

            if (pp.nPoints() && procPatch.owner())
            {
                // Get data per patchPoint in neighbouring point numbers.
                pointField patchInfo(procPatch.nPoints(), nullValue);

                const labelList& meshPts = procPatch.meshPoints();
                const labelList& nbrPts = procPatch.neighbPoints();

                forAll(nbrPts, pointi)
                {
                    label nbrPointi = nbrPts[pointi];
                    if (nbrPointi >= 0 && nbrPointi < patchInfo.size())
                    {
                        patchInfo[nbrPointi] = points[meshPts[pointi]];
                    }
                }

                OPstream toNbr
                (
                    UPstream::commsTypes::buffered,
                    procPatch.neighbProcNo()
                );
                toNbr << patchInfo;
            }
        }


        // Receive and set.

        for (const label patchi : procPatches)
        {
            const polyPatch& pp = patches[patchi];
            const auto& procPatch = refCast<const processorPolyPatch>(pp);

            if (pp.nPoints() && !procPatch.owner())
            {
                // We do not know the number of points on the other side
                // so cannot use UIPstream::read

                pointField nbrPatchInfo;
                IPstream::recv(nbrPatchInfo, procPatch.neighbProcNo());

                // Null any value which is not on neighbouring processor
                nbrPatchInfo.setSize(procPatch.nPoints(), nullValue);

                if (!procPatch.parallel())
                {
                    hasTransformation = true;
                    transformList(procPatch.forwardT(), nbrPatchInfo);
                }
                else if (procPatch.separated())
                {
                    hasTransformation = true;
                    separateList(-procPatch.separation(), nbrPatchInfo);
                }

                const labelList& meshPts = procPatch.meshPoints();

                forAll(meshPts, pointi)
                {
                    label meshPointi = meshPts[pointi];
                    points[meshPointi] = nbrPatchInfo[pointi];
                }
            }
        }
    }

    // Do the cyclics.
    for (const polyPatch& pp : patches)
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            // Owner does all.

            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            const edgeList& coupledPoints = cycPatch.coupledPoints();
            const labelList& meshPts = cycPatch.meshPoints();
            const labelList& nbrMeshPts = nbrPatch.meshPoints();

            pointField half0Values(coupledPoints.size());

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point0 = meshPts[e[0]];
                half0Values[i] = points[point0];
            }

            if (!cycPatch.parallel())
            {
                hasTransformation = true;
                transformList(cycPatch.reverseT(), half0Values);
            }
            else if (cycPatch.separated())
            {
                hasTransformation = true;
                separateList(cycPatch.separation(), half0Values);
            }

            forAll(coupledPoints, i)
            {
                const edge& e = coupledPoints[i];
                label point1 = nbrMeshPts[e[1]];
                points[point1] = half0Values[i];
            }
        }
    }

    //- Note: hasTransformation is only used for warning messages so
    //  reduction not strictly necessary.
    //Pstream::reduceOr(hasTransformation);

    // Synchronize multiple shared points.
    const globalMeshData& pd = mesh.globalData();

    if (pd.nGlobalPoints() > 0)
    {
        if (hasTransformation)
        {
            WarningInFunction
                << "There are decomposed cyclics in this mesh with"
                << " transformations." << endl
                << "This is not supported. The result will be incorrect"
                << endl;
        }


        // Values on shared points.
        pointField sharedPts(pd.nGlobalPoints(), nullValue);

        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            // Fill my entries in the shared points
            sharedPts[pd.sharedPointAddr()[i]] = points[meshPointi];
        }

        // Combine - globally consistent
        Pstream::listCombineReduce(sharedPts, cop);

        // Now we will all have the same information. Merge it back with
        // my local information.
        forAll(pd.sharedPointLabels(), i)
        {
            label meshPointi = pd.sharedPointLabels()[i];
            points[meshPointi] = sharedPts[pd.sharedPointAddr()[i]];
        }
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create patches out of selected boundary faces, which are either"
        " from existing patches or from a faceSet"
    );

    #include "addOverwriteOption.H"
    #include "addAllRegionOptions.H"

    argList::addOption("dict", "file", "Alternative createPatchDict");
    argList::addBoolOption
    (
        "writeObj",
        "Write obj files showing the cyclic matching process"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "getAllRegionOptions.H"

    const bool overwrite = args.found("overwrite");

    #include "createNamedMeshes.H"

    const bool writeObj = args.found("writeObj");


    // Read dictionaries and extract various
    wordList oldInstances(meshes.size());
    PtrList<dictionary> dicts(meshes.size());
    List<wordList> patchNames(meshes.size());
    List<labelListList> matchPatchIDs(meshes.size());
    List<wordList> matchMethods(meshes.size());
    List<PtrList<dictionary>> patchInfoDicts(meshes.size());

    List<DynamicList<label>> interRegionSources(meshes.size());
    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // If running parallel check same patches everywhere
        patches.checkParallelSync(true);

        oldInstances[meshi] = mesh.pointsInstance();

        const word dictName("createPatchDict");
        #include "setSystemMeshDictionaryIO.H"
        Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;

        dicts.set(meshi, new IOdictionary(dictIO));

        const auto& dict = dicts[meshi];
        PtrList<dictionary> patchSources(dict.lookup("patches"));
        patchNames[meshi].setSize(patchSources.size());
        matchPatchIDs[meshi].setSize(patchSources.size());
        matchMethods[meshi].setSize(patchSources.size());
        patchInfoDicts[meshi].setSize(patchSources.size());

        // Read patch construct info from dictionary
        forAll(patchSources, sourcei)
        {
            const auto& pDict = patchSources[sourcei];
            patchNames[meshi][sourcei] = pDict.getOrDefault<word>
            (
                "name",
                word::null,
                keyType::LITERAL
            );

            patchInfoDicts[meshi].set
            (
                sourcei,
                new dictionary(pDict.subDict("patchInfo"))
            );
            dictionary& patchDict = patchInfoDicts[meshi][sourcei];
            if (patchDict.found("AMIMethod"))
            {
                matchMethods[meshi][sourcei] = patchDict.get<word>("AMIMethod");
                // Disable full matching since we're trying to use AMIMethod to
                // find out actual overlap
                patchDict.add("requireMatch", false);
            }

            wordRes matchNames;
            if (pDict.readIfPresent("patches", matchNames, keyType::LITERAL))
            {
                matchPatchIDs[meshi][sourcei] =
                    patches.patchSet(matchNames).sortedToc();
            }

            if (pDict.get<word>("constructFrom") == "autoPatch")
            {
                interRegionSources[meshi].push_back(sourcei);
            }
        }
    }



    // Determine matching faces. This is an interface between two
    // sets of faces:
    //  - originating mesh
    //  - originating patch
    //  - originating (mesh)facelabels
    // It matches all mesh against each other. Lower numbered mesh gets
    // postfix 0, higher numbered mesh postfix 1.

    // Per interface, per patchSource:
    //      1. the lower numbered mesh
    DynamicList<labelList> interfaceMesh0;
    //      1b. the source index (i.e. the patch dictionary)
    DynamicList<label> interfaceSource0;
    //      2. the patch on the interfaceMesh0
    DynamicList<labelList> interfacePatch0;
    //      3. the facelabels on the interfaceMesh0
    DynamicList<List<DynamicList<label>>> interfaceFaces0;
    //      4. generated interface name
    DynamicList<wordList> interfaceNames0;

    // Same for the higher numbered mesh
    DynamicList<labelList> interfaceMesh1;
    DynamicList<label> interfaceSource1;
    DynamicList<labelList> interfacePatch1;
    DynamicList<List<DynamicList<label>>> interfaceFaces1;
    DynamicList<wordList> interfaceNames1;
    {
        // Whether to match to patches in own mesh
        const bool includeOwn = (meshes.size() == 1);

        matchPatchFaces
        (
            includeOwn,
            meshes,
            interRegionSources,
            patchNames,
            matchPatchIDs,
            matchMethods,
            patchInfoDicts,

            interfaceMesh0,
            interfaceSource0,
            interfacePatch0,
            interfaceFaces0,
            interfaceNames0,

            interfaceMesh1,
            interfaceSource1,
            interfacePatch1,
            interfaceFaces1,
            interfaceNames1
        );
    }



    // Read fields
    List<PtrList<volScalarField>> vsFlds(meshes.size());
    List<PtrList<volVectorField>> vvFlds(meshes.size());
    List<PtrList<volSphericalTensorField>> vstFlds(meshes.size());
    List<PtrList<volSymmTensorField>> vsymtFlds(meshes.size());
    List<PtrList<surfaceScalarField>> ssFlds(meshes.size());
    List<PtrList<volTensorField>> vtFlds(meshes.size());
    List<PtrList<surfaceVectorField>> svFlds(meshes.size());
    List<PtrList<surfaceSphericalTensorField>> sstFlds(meshes.size());
    List<PtrList<surfaceSymmTensorField>> ssymtFlds(meshes.size());
    List<PtrList<surfaceTensorField>> stFlds(meshes.size());

    {
        forAll(meshes, meshi)
        {
            const fvMesh& mesh = meshes[meshi];

            bool noFields = true;
            for (const auto& d : patchInfoDicts[meshi])
            {
                if (d.found("patchFields"))
                {
                    noFields = false;
                }
            }

            if (!noFields)
            {
                // Read objects in time directory
                IOobjectList objects(mesh, runTime.timeName());

                // Read volume fields.
                ReadFields(mesh, objects, vsFlds[meshi]);
                ReadFields(mesh, objects, vvFlds[meshi]);
                ReadFields(mesh, objects, vstFlds[meshi]);
                ReadFields(mesh, objects, vsymtFlds[meshi]);
                ReadFields(mesh, objects, vtFlds[meshi]);

                // Read surface fields.
                ReadFields(mesh, objects, ssFlds[meshi]);
                ReadFields(mesh, objects, svFlds[meshi]);
                ReadFields(mesh, objects, sstFlds[meshi]);
                ReadFields(mesh, objects, ssymtFlds[meshi]);
                ReadFields(mesh, objects, stFlds[meshi]);
            }
        }
    }



    // Loop over all regions

    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];

        Info<< "\n\nAdding patches to mesh " << mesh.name() << nl << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const dictionary& dict = dicts[meshi];

        if (writeObj)
        {
            dumpCyclicMatch("initial_", mesh);
        }

        // Read patch construct info from dictionary
        PtrList<dictionary> patchSources(dict.lookup("patches"));


        // 1. Add all new patches
        // ~~~~~~~~~~~~~~~~~~~~~~

        forAll(patchSources, sourcei)
        {
            const dictionary& dict = patchSources[sourcei];
            const word sourceType(dict.get<word>("constructFrom"));
            const word patchName
            (
                dict.getOrDefault<word>
                (
                    "name",
                    word::null,
                    keyType::LITERAL
                )
            );

            dictionary patchDict(patchInfoDicts[meshi][sourcei]);
            patchDict.set("nFaces", 0);
            patchDict.set("startFace", 0);    // Gets overwritten


            if (sourceType == "autoPatch")
            {
                // Special : automatically create the necessary inter-region
                // coupling patches

                forAll(interfaceMesh0, inti)
                {
                    const labelList& allMeshes0 = interfaceMesh0[inti];
                    const wordList& allNames0 = interfaceNames0[inti];
                    const labelList& allMeshes1 = interfaceMesh1[inti];
                    const wordList& allNames1 = interfaceNames1[inti];

                    if
                    (
                        interfaceSource0[inti] == sourcei
                     && allMeshes0[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh0. mesh1 is the remote mesh.
                        const label sourcej = interfaceSource1[inti];
                        const word& patchName = allNames0[sourcei];
                        if (patches.findPatchID(patchName) == -1)
                        {
                            dictionary allDict(patchDict);
                            const auto& mesh1 = meshes[allMeshes1[sourcej]];
                            allDict.set("sampleRegion", mesh1.name());
                            const auto& destPatch = allNames1[sourcej];
                            allDict.set("samplePatch", destPatch);
                            allDict.set("neighbourPatch", destPatch);

                            Info<< "Adding new patch " << patchName
                                << " from " << allDict << endl;

                            autoPtr<polyPatch> ppPtr
                            (
                                polyPatch::New
                                (
                                    patchName,
                                    allDict,
                                    0,          // overwritten
                                    patches
                                )
                            );
                            fvMeshTools::addPatch
                            (
                                mesh,
                                ppPtr(),
                                patchDict.subOrEmptyDict("patchFields"),
                                fvPatchFieldBase::calculatedType(),
                                true
                            );
                        }
                    }
                }

                forAll(interfaceMesh1, inti)
                {
                    const labelList& allMeshes0 = interfaceMesh0[inti];
                    const wordList& allNames0 = interfaceNames0[inti];
                    const labelList& allMeshes1 = interfaceMesh1[inti];
                    const wordList& allNames1 = interfaceNames1[inti];

                    if
                    (
                        interfaceSource1[inti] == sourcei
                     && allMeshes1[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh1. mesh0 is the remote mesh.

                        const label sourcej = interfaceSource0[inti];
                        const word& patchName = allNames1[sourcei];
                        if (patches.findPatchID(patchName) == -1)
                        {
                            dictionary allDict(patchDict);
                            const auto& destPatch = allNames0[sourcej];
                            const auto& mesh0 = meshes[allMeshes0[sourcej]];
                            allDict.set("sampleRegion", mesh0.name());
                            allDict.set("samplePatch", destPatch);
                            allDict.set("neighbourPatch", destPatch);

                            Info<< "Adding new patch " << patchName
                                << " from " << allDict << endl;

                            autoPtr<polyPatch> ppPtr
                            (
                                polyPatch::New
                                (
                                    patchName,
                                    allDict,
                                    0,          // overwritten
                                    patches
                                )
                            );
                            fvMeshTools::addPatch
                            (
                                mesh,
                                ppPtr(),
                                patchDict.subOrEmptyDict("patchFields"),
                                fvPatchFieldBase::calculatedType(),
                                true
                            );
                        }
                    }
                }
            }
            else
            {
                if (patches.findPatchID(patchName) == -1)
                {
                    Info<< "Adding new patch " << patchName
                        << " from " << patchDict << endl;

                    autoPtr<polyPatch> ppPtr
                    (
                        polyPatch::New
                        (
                            patchName,
                            patchDict,
                            0,          // overwritten
                            patches
                        )
                    );
                    fvMeshTools::addPatch
                    (
                        mesh,
                        ppPtr(),
                        patchDict.subOrEmptyDict("patchFields"),
                        fvPatchFieldBase::calculatedType(),
                        true
                    );
                }
            }
        }

        Info<< endl;
    }


    // 2. Repatch faces
    // ~~~~~~~~~~~~~~~~

    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];

        Info<< "\n\nRepatching mesh " << mesh.name() << nl << endl;


        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const dictionary& dict = dicts[meshi];

        // Whether to synchronise points
        const bool pointSync(dict.get<bool>("pointSync"));

        // Read patch construct info from dictionary
        PtrList<dictionary> patchSources(dict.lookup("patches"));


        polyTopoChange meshMod(mesh);

        // Mark all repatched faces. This makes sure that the faces to repatch
        // do not overlap
        bitSet isRepatchedBoundary(mesh.nBoundaryFaces());

        forAll(patchSources, sourcei)
        {
            const dictionary& dict = patchSources[sourcei];
            const word patchName
            (
                dict.getOrDefault<word>
                (
                    "name",
                    word::null,
                    keyType::LITERAL
                )
            );
            const word sourceType(dict.get<word>("constructFrom"));

            if (sourceType == "autoPatch")
            {
                forAll(interfaceMesh0, inti)
                {
                    const labelList& allMeshes0 = interfaceMesh0[inti];
                    const wordList& allNames0 = interfaceNames0[inti];
                    const auto& allFaces0 = interfaceFaces0[inti];

                    const labelList& allMeshes1 = interfaceMesh1[inti];
                    const wordList& allNames1 = interfaceNames1[inti];
                    const auto& allFaces1 = interfaceFaces1[inti];

                    if
                    (
                        interfaceSource0[inti] == sourcei
                     && allMeshes0[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh0. mesh1 is the remote mesh.

                        const label destPatchi =
                            patches.findPatchID(allNames0[sourcei], false);

                        //const auto& mesh1 =
                        //    meshes[allMeshes1[interfaceSource1[inti]]];
                        //Pout<< "Matched mesh:" << mesh.name()
                        //    << " to mesh:" << mesh1.name()
                        //    << " through:" << allNames0[sourcei] << endl;

                        changePatchID
                        (
                            mesh,
                            allFaces0[sourcei],
                            destPatchi,
                            isRepatchedBoundary,
                            meshMod
                        );
                    }
                    if
                    (
                        interfaceSource1[inti] == sourcei
                     && allMeshes1[sourcei] == meshi
                    )
                    {
                        // Current mesh is mesh1. mesh0 is the remote mesh.

                        const label destPatchi =
                            patches.findPatchID(allNames1[sourcei], false);

                        //const auto& mesh0 =
                        //    meshes[allMeshes0[interfaceSource0[inti]]];
                        //Pout<< "Matched mesh:" << mesh.name()
                        //    << " to mesh:" << mesh0.name()
                        //    << " through:" << allNames1[sourcei] << endl;

                        changePatchID
                        (
                            mesh,
                            allFaces1[sourcei],
                            destPatchi,
                            isRepatchedBoundary,
                            meshMod
                        );
                    }
                }
            }
            else if (sourceType == "patches")
            {
                const label destPatchi = patches.findPatchID(patchName, false);
                labelHashSet patchSources
                (
                    patches.patchSet(dict.get<wordRes>("patches"))
                );

                // Repatch faces of the patches.
                for (const label patchi : patchSources)
                {
                    const polyPatch& pp = patches[patchi];

                    Info<< "Moving faces from patch " << pp.name()
                        << " to patch " << destPatchi << endl;

                    changePatchID
                    (
                        mesh,
                        pp.start()+identity(pp.size()),
                        destPatchi,
                        isRepatchedBoundary,
                        meshMod
                    );
                }
            }
            else if (sourceType == "set")
            {
                const label destPatchi = patches.findPatchID(patchName, false);
                const word setName(dict.get<word>("set"));

                faceSet set(mesh, setName);

                Info<< "Read " << returnReduce(set.size(), sumOp<label>())
                    << " faces from faceSet " << set.name() << endl;

                // Sort (since faceSet contains faces in arbitrary order)
                labelList faceLabels(set.sortedToc());

                changePatchID
                (
                    mesh,
                    faceLabels,
                    destPatchi,
                    isRepatchedBoundary,
                    meshMod
                );
            }
            else
            {
                FatalErrorInFunction
                    << "Invalid source type " << sourceType << endl
                    << "Valid source types are 'patches' 'set'"
                    << exit(FatalError);
            }
        }


        // Change mesh, use inflation to reforce calculation of transformation
        // tensors.
        Info<< "Doing topology modification to order faces." << nl << endl;
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);

        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
           // Force calculation of transformation tensors
            mesh.movePoints(pointField(mesh.points()));
        }


        // Update fields
        mesh.updateMesh(map());


        // Update numbering pointing to meshi
        const auto& oldToNew = map().reverseFaceMap();
        forAll(interfaceMesh0, inti)
        {
            const labelList& allMeshes0 = interfaceMesh0[inti];
            forAll(allMeshes0, sourcei)
            {
                if (allMeshes0[sourcei] == meshi)
                {
                    inplaceRenumber(oldToNew, interfaceFaces0[inti][sourcei]);
                }
            }
        }
        forAll(interfaceMesh1, inti)
        {
            const labelList& allMeshes1 = interfaceMesh1[inti];
            forAll(allMeshes1, sourcei)
            {
                if (allMeshes1[sourcei] == meshi)
                {
                    inplaceRenumber(oldToNew, interfaceFaces1[inti][sourcei]);
                }
            }
        }


        if (writeObj)
        {
            dumpCyclicMatch("coupled_", mesh);
        }

        // Synchronise points.
        if (!pointSync)
        {
            Info<< "Not synchronising points." << nl << endl;
        }
        else
        {
            Info<< "Synchronising points." << nl << endl;

            // This is a bit tricky. Both normal and position might be out and
            // current separation also includes the normal
            // ( separation_ = (nf&(Cr - Cf))*nf ).

            // For cyclic patches:
            // - for separated ones use user specified offset vector

            forAll(mesh.boundaryMesh(), patchi)
            {
                const polyPatch& pp = mesh.boundaryMesh()[patchi];

                if (pp.size() && isA<coupledPolyPatch>(pp))
                {
                    const coupledPolyPatch& cpp =
                        refCast<const coupledPolyPatch>(pp);

                    if (cpp.separated())
                    {
                        Info<< "On coupled patch " << pp.name()
                            << " separation[0] was "
                            << cpp.separation()[0] << endl;

                        if (isA<cyclicPolyPatch>(pp) && pp.size())
                        {
                            const auto& cycpp =
                                refCast<const cyclicPolyPatch>(pp);

                            if
                            (
                                cycpp.transform()
                             == cyclicPolyPatch::TRANSLATIONAL
                            )
                            {
                                // Force to wanted separation
                                Info<< "On cyclic translation patch "
                                    << pp.name()
                                    << " forcing uniform separation of "
                                    << cycpp.separationVector() << endl;
                                const_cast<vectorField&>(cpp.separation()) =
                                    pointField(1, cycpp.separationVector());
                            }
                            else
                            {
                                const auto& nbr = cycpp.neighbPatch();
                                const_cast<vectorField&>(cpp.separation()) =
                                    pointField
                                    (
                                        1,
                                        nbr[0].centre(mesh.points())
                                      - cycpp[0].centre(mesh.points())
                                    );
                            }
                        }
                        Info<< "On coupled patch " << pp.name()
                            << " forcing uniform separation of "
                            << cpp.separation() << endl;
                    }
                    else if (!cpp.parallel())
                    {
                        Info<< "On coupled patch " << pp.name()
                            << " forcing uniform rotation of "
                            << cpp.forwardT()[0] << endl;

                        const_cast<tensorField&>
                        (
                            cpp.forwardT()
                        ).setSize(1);
                        const_cast<tensorField&>
                        (
                            cpp.reverseT()
                        ).setSize(1);

                        Info<< "On coupled patch " << pp.name()
                            << " forcing uniform rotation of "
                            << cpp.forwardT() << endl;
                    }
                }
            }

            Info<< "Synchronising points." << endl;

            pointField newPoints(mesh.points());

            syncPoints
            (
                mesh,
                newPoints,
                minMagSqrEqOp<vector>(),
                point(GREAT, GREAT, GREAT)
            );

            scalarField diff(mag(newPoints-mesh.points()));
            Info<< "Points changed by average:" << gAverage(diff)
                << " max:" << gMax(diff) << nl << endl;

            mesh.movePoints(newPoints);
        }


        // 3. Remove zeros-sized patches
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Info<< "Removing patches with no faces in them." << nl << endl;
        const wordList oldPatchNames(mesh.boundaryMesh().names());
        const wordList oldPatchTypes(mesh.boundaryMesh().types());
        fvMeshTools::removeEmptyPatches(mesh, true);
        forAll(oldPatchNames, patchi)
        {
            const word& pName = oldPatchNames[patchi];
            if (mesh.boundaryMesh().findPatchID(pName) == -1)
            {
                Info<< "Removed zero-sized patch " << pName
                    << " type " << oldPatchTypes[patchi]
                    << " at position " << patchi << endl;
            }
        }


        if (writeObj)
        {
            dumpCyclicMatch("final_", mesh);
        }
    }

    if (!overwrite)
    {
        ++runTime;
    }
    else
    {
        forAll(meshes, meshi)
        {
            fvMesh& mesh = meshes[meshi];

            mesh.setInstance(oldInstances[meshi]);
        }
    }

    // More precision (for points data)
    IOstream::minPrecision(10);

    // Write resulting mesh
    forAll(meshes, meshi)
    {
        fvMesh& mesh = meshes[meshi];
        Info<< "\n\nWriting repatched mesh " << mesh.name()
            << " to " << runTime.timeName() << nl << endl;
        mesh.clearOut();    // remove meshPhi
        mesh.write();
        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

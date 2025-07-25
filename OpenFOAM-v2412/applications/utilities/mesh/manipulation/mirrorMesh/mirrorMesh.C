/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    mirrorMesh

Group
    grpMeshManipulationUtilities

Description
    Mirrors a mesh around a given plane.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "mirrorFvMesh.H"
#include "mapPolyMesh.H"
#include "hexRef8Data.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mirrors a mesh around a given plane."
    );

    argList::addOption("dict", "file", "Alternative mirrorMeshDict");
    argList::setAdvanced("decomposeParDict");

    #include "addOverwriteOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const bool overwrite = args.found("overwrite");

    const word dictName("mirrorMeshDict");

    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictIO.name() << nl << endl;

    const IOdictionary mirrorDict(dictIO);

    mirrorFvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime
        ),
        mirrorDict
    );

    hexRef8Data refData
    (
        IOobject
        (
            "dummy",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    if (!overwrite)
    {
        ++runTime;
        mesh.setInstance(runTime.timeName());
    }


    // More precision (for points data)
    IOstream::minPrecision(10);

    // Generate the mirrorred mesh
    const fvMesh& mMesh = mesh.mirrorMesh();

    const_cast<fvMesh&>(mMesh).setInstance(mesh.facesInstance());
    Info<< "Writing mirrored mesh" << endl;
    mMesh.write();

    // Map the hexRef8 data
    mapPolyMesh map
    (
        mesh,
        mesh.nPoints(),         //nOldPoints,
        mesh.nFaces(),          //nOldFaces,
        mesh.nCells(),          //nOldCells,
        mesh.pointMap(),        //pointMap,
        List<objectMap>(0),     // pointsFromPoints,
        labelList(0),           //faceMap,
        List<objectMap>(0),     //facesFromPoints,
        List<objectMap>(0),     //facesFromEdges,
        List<objectMap>(0),     //facesFromFaces,
        mesh.cellMap(),         //cellMap,
        List<objectMap>(0),     //cellsFromPoints,
        List<objectMap>(0),     //cellsFromEdges,
        List<objectMap>(0),     //cellsFromFaces,
        List<objectMap>(0),     //cellsFromCells,
        labelList(0),           //reversePointMap,
        labelList(0),           //reverseFaceMap,
        labelList(0),           //reverseCellMap,
        labelHashSet(0),        //flipFaceFlux,
        labelListList(0),       //patchPointMap,
        labelListList(0),       //pointZoneMap,
        labelListList(0),       //faceZonePointMap,
        labelListList(0),       //faceZoneFaceMap,
        labelListList(0),       //cellZoneMap,
        pointField(0),          //preMotionPoints,
        labelList(0),           //oldPatchStarts,
        labelList(0),           //oldPatchNMeshPoints,
        autoPtr<scalarField>()  //oldCellVolumesPtr
    );
    refData.updateMesh(map);
    refData.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

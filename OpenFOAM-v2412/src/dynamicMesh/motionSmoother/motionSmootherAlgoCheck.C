/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2020,2022 OpenCFD Ltd.
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

#include "motionSmootherAlgo.H"
#include "polyMeshGeometry.H"
#include "IOmanip.H"
#include "pointSet.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    labelHashSet& wrongFaces,
    const bool dryRun
)
{
    List<labelPair> emptyBaffles;
    return checkMesh
    (
        report,
        mesh,
        dict,
        checkFaces,
        emptyBaffles,
        wrongFaces,
        dryRun
    );
}

bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces,
    const bool dryRun
)
{
    const scalar maxNonOrtho
    (
        get<scalar>(dict, "maxNonOrtho", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minVol
    (
        get<scalar>(dict, "minVol", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minTetQuality
    (
        get<scalar>(dict, "minTetQuality", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar maxConcave
    (
        get<scalar>(dict, "maxConcave", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minArea
    (
        get<scalar>(dict, "minArea", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar maxIntSkew
    (
        get<scalar>
        (
            dict, "maxInternalSkewness", dryRun, keyType::REGEX_RECURSIVE
        )
    );
    const scalar maxBounSkew
    (
        get<scalar>
        (
            dict, "maxBoundarySkewness", dryRun, keyType::REGEX_RECURSIVE
        )
    );
    const scalar minWeight
    (
        get<scalar>(dict, "minFaceWeight", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minVolRatio
    (
        get<scalar>(dict, "minVolRatio", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minTwist
    (
        get<scalar>(dict, "minTwist", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minTriangleTwist
    (
        get<scalar>(dict, "minTriangleTwist", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minFaceFlatness
    (
        dict.getOrDefault<scalar>
        (
            "minFaceFlatness", -1, keyType::REGEX_RECURSIVE
        )
    );
    const scalar minDet
    (
        get<scalar>(dict, "minDeterminant", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minEdgeLength
    (
        dict.getOrDefault<scalar>
        (
            "minEdgeLength", -1, keyType::REGEX_RECURSIVE
        )
    );


    if (dryRun)
    {
        string errorMsg(FatalError.message());
        string IOerrorMsg(FatalIOError.message());

        if (errorMsg.size() || IOerrorMsg.size())
        {
            //errorMsg = "[dryRun] " + errorMsg;
            //errorMsg.replaceAll("\n", "\n[dryRun] ");
            //IOerrorMsg = "[dryRun] " + IOerrorMsg;
            //IOerrorMsg.replaceAll("\n", "\n[dryRun] ");

            IOWarningInFunction(dict)
                << nl
                << "Missing/incorrect required dictionary entries:" << nl
                << nl
                << IOerrorMsg.c_str() << nl
                << errorMsg.c_str() << nl
                //<< nl << "Exiting dry-run" << nl
                << endl;

            FatalError.clear();
            FatalIOError.clear();
        }
        return false;
    }


    label nWrongFaces = 0;

    Info<< "Checking faces in error :" << endl;
    //Pout.setf(ios_base::left);

    if (maxNonOrtho < 180.0-SMALL)
    {
        polyMeshGeometry::checkFaceDotProduct
        (
            report,
            maxNonOrtho,
            mesh,
            mesh.cellCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    non-orthogonality > "
            << setw(3) << maxNonOrtho
            << " degrees                        : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVol > -GREAT)
    {
        polyMeshGeometry::checkFacePyramids
        (
            report,
            minVol,
            mesh,
            mesh.cellCentres(),
            mesh.points(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face pyramid volume < "
            << setw(5) << minVol << "                 : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTetQuality > -GREAT)
    {
        polyMeshGeometry::checkFaceTets
        (
            report,
            minTetQuality,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face-decomposition tet quality < "
            << setw(5) << minTetQuality << "      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxConcave < 180.0-SMALL)
    {
        polyMeshGeometry::checkFaceAngles
        (
            report,
            maxConcave,
            mesh,
            mesh.faceAreas(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with concavity > "
            << setw(3) << maxConcave
            << " degrees                     : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minArea > -SMALL)
    {
        polyMeshGeometry::checkFaceArea
        (
            report,
            minArea,
            mesh,
            mesh.faceAreas(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with area < "
            << setw(5) << minArea
            << " m^2                            : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxIntSkew > 0 || maxBounSkew > 0)
    {
        polyMeshGeometry::checkFaceSkewness
        (
            report,
            maxIntSkew,
            maxBounSkew,
            mesh,
            mesh.points(),
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with skewness > "
            << setw(3) << maxIntSkew
            << " (internal) or " << setw(3) << maxBounSkew
            << " (boundary) : " << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minWeight >= 0 && minWeight < 1)
    {
        polyMeshGeometry::checkFaceWeights
        (
            report,
            minWeight,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with interpolation weights (0..1)  < "
            << setw(5) << minWeight
            << "       : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVolRatio >= 0)
    {
        polyMeshGeometry::checkVolRatio
        (
            report,
            minVolRatio,
            mesh,
            mesh.cellVolumes(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with volume ratio of neighbour cells < "
            << setw(5) << minVolRatio
            << "     : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTwist > -1)
    {
        //Pout<< "Checking face twist: dot product of face normal "
        //    << "with face triangle normals" << endl;
        polyMeshGeometry::checkFaceTwist
        (
            report,
            minTwist,
            mesh,
            mesh.cellCentres(),
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face twist < "
            << setw(5) << minTwist
            << "                          : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTriangleTwist > -1)
    {
        //Pout<< "Checking triangle twist: dot product of consecutive triangle"
        //    << " normals resulting from face-centre decomposition" << endl;
        polyMeshGeometry::checkTriangleTwist
        (
            report,
            minTriangleTwist,
            mesh,
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with triangle twist < "
            << setw(5) << minTriangleTwist
            << "                      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minFaceFlatness > -SMALL)
    {
        polyMeshGeometry::checkFaceFlatness
        (
            report,
            minFaceFlatness,
            mesh,
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with flatness < "
            << setw(5) << minFaceFlatness
            << "                      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minDet > -1)
    {
        polyMeshGeometry::checkCellDeterminant
        (
            report,
            minDet,
            mesh,
            mesh.faceAreas(),
            checkFaces,
            polyMeshGeometry::affectedCells(mesh, checkFaces),
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces on cells with determinant < "
            << setw(5) << minDet << "                : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minEdgeLength >= 0)
    {
        pointSet edgePoints(mesh, "smallEdgePoints", mesh.nPoints()/1000);
        polyMeshGeometry::checkEdgeLength
        (
            report,
            minEdgeLength,
            mesh,
            checkFaces,
            &edgePoints
        );

        const auto& pointFaces = mesh.pointFaces();

        label nNewWrongFaces = 0;
        for (const label pointi : edgePoints)
        {
            const auto& pFaces = pointFaces[pointi];
            for (const label facei : pFaces)
            {
                if (wrongFaces.insert(facei))
                {
                    nNewWrongFaces++;
                }
            }
        }

        Info<< "    faces with edge length < "
            << setw(5) << minEdgeLength << "                         : "
            << returnReduce(nNewWrongFaces, sumOp<label>()) << endl;

        nWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());
    }

    //Pout.setf(ios_base::right);

    return nWrongFaces > 0;
}


bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    labelHashSet& wrongFaces,
    const bool dryRun
)
{
    return checkMesh
    (
        report,
        mesh,
        dict,
        identity(mesh.nFaces()),
        wrongFaces,
        dryRun
    );
}

bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const dictionary& dict,
    const polyMeshGeometry& meshGeom,
    const pointField& points,
    const labelList& checkFaces,
    labelHashSet& wrongFaces,
    const bool dryRun
)
{
    List<labelPair> emptyBaffles;

    return checkMesh
    (
        report,
        dict,
        meshGeom,
        points,
        checkFaces,
        emptyBaffles,
        wrongFaces,
        dryRun
     );
}


bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const dictionary& dict,
    const polyMeshGeometry& meshGeom,
    const pointField& points,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces,
    const bool dryRun
)
{
    const scalar maxNonOrtho
    (
        get<scalar>(dict, "maxNonOrtho", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minVol
    (
        get<scalar>(dict, "minVol", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minTetQuality
    (
        get<scalar>(dict, "minTetQuality", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar maxConcave
    (
        get<scalar>(dict, "maxConcave", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minArea
    (
        get<scalar>(dict, "minArea", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar maxIntSkew
    (
        get<scalar>
        (
            dict,
            "maxInternalSkewness",
            dryRun,
            keyType::REGEX_RECURSIVE
        )
    );
    const scalar maxBounSkew
    (
        get<scalar>
        (
            dict,
            "maxBoundarySkewness",
            dryRun,
            keyType::REGEX_RECURSIVE
        )
    );
    const scalar minWeight
    (
        get<scalar>(dict, "minFaceWeight", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minVolRatio
    (
        get<scalar>(dict, "minVolRatio", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minTwist
    (
        get<scalar>(dict, "minTwist", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minTriangleTwist
    (
        get<scalar>(dict, "minTriangleTwist", dryRun, keyType::REGEX_RECURSIVE)
    );
    scalar minFaceFlatness = -1.0;
    dict.readIfPresent
    (
        "minFaceFlatness",
        minFaceFlatness,
        keyType::REGEX_RECURSIVE
    );
    const scalar minDet
    (
        get<scalar>(dict, "minDeterminant", dryRun, keyType::REGEX_RECURSIVE)
    );
    const scalar minEdgeLength
    (
        dict.getOrDefault<scalar>
        (
            "minEdgeLength", -1, keyType::REGEX_RECURSIVE
        )
    );

    if (dryRun)
    {
        string errorMsg(FatalError.message());
        string IOerrorMsg(FatalIOError.message());

        if (errorMsg.size() || IOerrorMsg.size())
        {
            //errorMsg = "[dryRun] " + errorMsg;
            //errorMsg.replaceAll("\n", "\n[dryRun] ");
            //IOerrorMsg = "[dryRun] " + IOerrorMsg;
            //IOerrorMsg.replaceAll("\n", "\n[dryRun] ");

            Perr<< nl
                << "Missing/incorrect required dictionary entries:" << nl
                << nl
                << IOerrorMsg.c_str() << nl
                << errorMsg.c_str() << nl
                //<< nl << "Exiting dry-run" << nl
                << endl;

            FatalError.clear();
            FatalIOError.clear();
        }
        return false;
    }


    label nWrongFaces = 0;

    Info<< "Checking faces in error :" << endl;
    //Pout.setf(ios_base::left);

    if (maxNonOrtho < 180.0-SMALL)
    {
        meshGeom.checkFaceDotProduct
        (
            report,
            maxNonOrtho,
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    non-orthogonality > "
            << setw(3) << maxNonOrtho
            << " degrees                        : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVol > -GREAT)
    {
        meshGeom.checkFacePyramids
        (
            report,
            minTetQuality,
            points,
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face pyramid volume < "
            << setw(5) << minVol << "                 : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTetQuality > -GREAT)
    {
        meshGeom.checkFaceTets
        (
            report,
            minTetQuality,
            points,
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face-decomposition tet quality < "
            << setw(5) << minTetQuality << "      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxConcave < 180.0-SMALL)
    {
        meshGeom.checkFaceAngles
        (
            report,
            maxConcave,
            points,
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with concavity > "
            << setw(3) << maxConcave
            << " degrees                     : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minArea > -SMALL)
    {
        meshGeom.checkFaceArea
        (
            report,
            minArea,
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with area < "
            << setw(5) << minArea
            << " m^2                            : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxIntSkew > 0 || maxBounSkew > 0)
    {
        polyMeshGeometry::checkFaceSkewness
        (
            report,
            maxIntSkew,
            maxBounSkew,
            meshGeom.mesh(),
            points,
            meshGeom.cellCentres(),
            meshGeom.faceCentres(),
            meshGeom.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with skewness > "
            << setw(3) << maxIntSkew
            << " (internal) or " << setw(3) << maxBounSkew
            << " (boundary) : " << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minWeight >= 0 && minWeight < 1)
    {
        meshGeom.checkFaceWeights
        (
            report,
            minWeight,
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with interpolation weights (0..1)  < "
            << setw(5) << minWeight
            << "       : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVolRatio >= 0)
    {
        meshGeom.checkVolRatio
        (
            report,
            minVolRatio,
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with volume ratio of neighbour cells < "
            << setw(5) << minVolRatio
            << "     : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTwist > -1)
    {
        //Pout<< "Checking face twist: dot product of face normal "
        //    << "with face triangle normals" << endl;
        meshGeom.checkFaceTwist
        (
            report,
            minTwist,
            points,
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face twist < "
            << setw(5) << minTwist
            << "                          : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTriangleTwist > -1)
    {
        //Pout<< "Checking triangle twist: dot product of consecutive triangle"
        //    << " normals resulting from face-centre decomposition" << endl;
        meshGeom.checkTriangleTwist
        (
            report,
            minTriangleTwist,
            points,
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with triangle twist < "
            << setw(5) << minTriangleTwist
            << "                      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minFaceFlatness > -SMALL)
    {
        meshGeom.checkFaceFlatness
        (
            report,
            minFaceFlatness,
            points,
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with flatness < "
            << setw(5) << minFaceFlatness
            << "                      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minDet > -1)
    {
        meshGeom.checkCellDeterminant
        (
            report,
            minDet,
            checkFaces,
            polyMeshGeometry::affectedCells(meshGeom.mesh(), checkFaces),
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces on cells with determinant < "
            << setw(5) << minDet << "                : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minEdgeLength >= 0)
    {
        pointSet edgePoints
        (
            meshGeom.mesh(),
            "smallEdgePoints",
            meshGeom.mesh().nPoints()/1000
        );
        meshGeom.checkEdgeLength
        (
            report,
            minEdgeLength,
            checkFaces,
            &edgePoints
        );

        const auto& pointFaces = meshGeom.mesh().pointFaces();

        label nNewWrongFaces = 0;
        for (const label pointi : edgePoints)
        {
            const auto& pFaces = pointFaces[pointi];
            for (const label facei : pFaces)
            {
                if (wrongFaces.insert(facei))
                {
                    nNewWrongFaces++;
                }
            }
        }

        Info<< "    faces with edge length < "
            << setw(5) << minEdgeLength << "                         : "
            << returnReduce(nNewWrongFaces, sumOp<label>()) << endl;

        nWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());
    }

    //Pout.setf(ios_base::right);

    return nWrongFaces > 0;
}


// ************************************************************************* //

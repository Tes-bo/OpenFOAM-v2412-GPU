/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021,2024 OpenCFD Ltd.
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

#include "hexMeshSmootherMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "pointConstraints.H"
#include "unitConversion.H"
#include "OBJstream.H"
#include "PatchTools.H"
//#include "geometricElementTransformPointSmoother.H"
#include "pointList.H"
#include "vectorList.H"
#include "meshPointPatch.H"
#include "pointSmoother.H"
#include "fvMesh.H"
#include "fvGeometryScheme.H"
#include "emptyPointPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(hexMeshSmootherMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        hexMeshSmootherMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        hexMeshSmootherMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::hexMeshSmootherMotionSolver::nonConstraintPatches
(
    const polyMesh& mesh
)
{
    // Get list of all non-constraint patches. These are the ones where
    // laplacian smoothing is applied.

    const auto& pbm = mesh.boundaryMesh();

    DynamicList<label> patchIDs(pbm.size());
    for (const auto& pp : pbm)
    {
        if (!polyPatch::constraintType(pp.type()))
        {
            patchIDs.append(pp.index());
        }
    }
    return patchIDs;
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::hexMeshSmootherMotionSolver::makePatch
(
    const polyMesh& mesh,
    const labelList& patchIDs,
    const labelList& zoneIDs,
    const pointField& points0
)
{
    // Mark all faces
    bitSet isPatchFace(mesh.nFaces());

    // Mark all boundary faces (or just patchIDs?)
    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        isPatchFace.set(pp.range());
    }

    const auto& fzs = mesh.faceZones();
    for (const label zonei : zoneIDs)
    {
        isPatchFace.set(fzs[zonei]);
    }

    syncTools::syncFaceList(mesh, isPatchFace, orEqOp<unsigned int>());

    const labelList patchFaces(isPatchFace.sortedToc());

    return autoPtr<indirectPrimitivePatch>::New
    (
        IndirectList<face>(mesh.faces(), patchFaces),
        points0
    );
}


void Foam::hexMeshSmootherMotionSolver::checkMesh
(
    const pointField& currentPoints,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs,
    const scalarField& cellVols,
    labelHashSet& markedFaces,
    bitSet& markedPoints
) const
{
    // Replacement for motionSmootherAlgo::checkMesh. Adds to markedFaces
    // any faces that are insufficient quality


    markedFaces.clear();
    markedPoints = false;


    /*

    tmp<scalarField> tminCellQ
    (
        pointSmoothers::geometricElementTransformPointSmoother::cellQuality
        (
            mesh(),
            currentPoints
        )
    );
    const scalarField& minCellQ = tminCellQ();

    markedPoints.setSize(mesh().nPoints());

    labelHashSet set;
    DynamicList<label> storage;

    for (label facei = 0; facei < mesh().nFaces(); facei++)
    {
        const label own = mesh().faceOwner()[facei];
        if (minCellQ[own] < VSMALL)
        {
            markedFaces.insert(facei);
            markedPoints.set(mesh().cellPoints(own, set, storage));
        }
        else if
        (
            mesh().isInternalFace(facei)
         && minCellQ[mesh().faceNeighbour()[facei]] < VSMALL
        )
        {
            markedFaces.insert(facei);
            markedPoints.set
            (
                mesh().cellPoints(mesh().faceNeighbour()[facei], set, storage)
            );
        }
    }
    */

    // Get measure of face quality
    tmp<scalarField> tfaceQ
    (
        pointSmoother_->faceQuality
        (
            currentPoints,
            fCtrs,
            fAreas,
            cellCtrs,
            cellVols
        )
    );
    const auto& faceQ = tfaceQ();

    markedPoints.setSize(mesh().nPoints());
    forAll(faceQ, facei)
    {
        if (faceQ[facei] < VSMALL)
        {
            markedFaces.insert(facei);
            markedPoints.set(mesh().faces()[facei]);
        }
    }


    syncTools::syncPointList
    (
        mesh(),
        markedPoints,
        orEqOp<unsigned int>(),
        0U
    );

    // Par sync. TBD.
    {
        bitSet isMarkedFace(mesh().nFaces());
        isMarkedFace.set(markedFaces.toc());
        syncTools::syncFaceList
        (
            mesh(),
            isMarkedFace,
            orEqOp<unsigned int>()
        );
        markedFaces.insert(isMarkedFace.toc());
    }
}


//void Foam::hexMeshSmootherMotionSolver::constrainDisplacement
//(
//    pointField& points
//) const
//{
//    // Make sure the points obey the boundary conditions
//    // on pointDisplacement
//
//    // Update pointDisplacement for suppled points
//    pointDisplacement_.primitiveFieldRef() = points-points0();
//    const pointConstraints& pcs =
//        pointConstraints::New(pointDisplacement_.mesh());
//    pcs.constrainDisplacement(pointDisplacement_, false);
////    pointDisplacement_.correctBoundaryConditions();
//    points = points0()+pointDisplacement();
//}


bool Foam::hexMeshSmootherMotionSolver::relax
(
    const scalarList& relaxationFactors,
    const bitSet& pointsToRelax,
    const pointField& initialPoints,
    const pointField& wantedPoints,
    pointField& relaxedPoints,
    labelList& relaxationLevel
) const
{
    // Find relaxation level that makes mesh quality acceptable. Gets given
    // initial set of mesh points

    relaxedPoints = wantedPoints;

    {
        vectorField fCtrs(mesh().nFaces());
        vectorField fAreas(mesh().nFaces());
        vectorField cellCtrs(mesh().nCells());
        scalarField cellVols(mesh().nCells());

        // Calculate mesh quantities with new locations
        const auto& geom =
            reinterpret_cast<const fvMesh&>(mesh()).geometry();
        geom.updateGeom
        (
            relaxedPoints,
            mesh().points(),    // old points (for avoiding recalculation)
            fCtrs,
            fAreas,
            cellCtrs,
            cellVols
        );

        // Check the modified face quality. Marks faces with insufficient
        // quality.
        labelHashSet markedFaces;
        bitSet markedPoints;
        checkMesh
        (
            relaxedPoints,
            fCtrs,
            fAreas,
            cellCtrs,
            cellVols,
            markedFaces,
            markedPoints
        );
        if (debug)
        {
            Pout<< "** hexMeshSmootherMotionSolver::relax : errorfaces:"
                << markedFaces.size()
                << " errorpoints:" << markedPoints.count() << endl;
        }
    }


    // Create a list of relaxation levels
    // -1 indicates a point which is not to be moved
    //  0 is the starting value for a moving point
    relaxationLevel.setSize(mesh().nPoints());
    relaxationLevel = -1;
    for (const label pointi : pointsToRelax)
    {
        relaxationLevel[pointi] = 0;
    }

    syncTools::syncPointList
    (
        mesh(),
        relaxationLevel,
        maxEqOp<label>(),
        label(-1)
    );

    vectorField fCtrs(mesh().nFaces());
    vectorField fAreas(mesh().nFaces());
    vectorField cellCtrs(mesh().nCells());
    scalarField cellVols(mesh().nCells());

    // Loop whilst relaxation levels are being incremented
    bool complete(false);
    while (!complete)
    {
        //Info<< "    Moving "
        //    << countZeroOrPos(relaxationFactors.size(), relaxationLevel)
        //    << " points" << endl;

        // Calculate current points (relaxationLevel >= 0)
        forAll(relaxationLevel, pointi)
        {
            if (relaxationLevel[pointi] >= 0)
            {
                const scalar x
                (
                    relaxationFactors[relaxationLevel[pointi]]
                );

                relaxedPoints[pointi] =
                    (1 - x)*initialPoints[pointi]
                  + x*wantedPoints[pointi];
            }
        }

        // Make sure the relaxed points still obey the boundary conditions
        // on pointDisplacement. Note: could do this afterwards but better
        // as soon as possible so we pick it up in the checkMesh
        //constrainDisplacement(relaxedPoints);


        // Calculate mesh quantities with new locations
        const auto& geom =
            reinterpret_cast<const fvMesh&>(mesh()).geometry();
        geom.updateGeom
        (
            relaxedPoints,
            mesh().points(),    // old points
            fCtrs,
            fAreas,
            cellCtrs,
            cellVols
        );

        // Check the modified face quality. Marks faces with insufficient
        // quality.
        labelHashSet markedFaces;
        bitSet markedPoints;
        checkMesh
        (
            relaxedPoints,
            fCtrs,
            fAreas,
            cellCtrs,
            cellVols,
            markedFaces,
            markedPoints
        );
        //Pout<< "    checkMesh : errorfaces:" << markedFaces.size()
        //    << " errorpoints:" << markedPoints.count() << endl;

        complete = true;
        for (const label pointi : markedPoints)
        {
            if (relaxationLevel[pointi] < relaxationFactors.size() - 1)
            {
                ++relaxationLevel[pointi];
                complete = false;
            }
        }

        //Info<< "    After adjustment:"
        //    << countZeroOrPos(relaxationFactors.size(), relaxationLevel)
        //    << " points relaxed" << endl;


        // Synchronise convergence
        reduce(complete, andOp<bool>());

        // Synchronise relaxation levels
        syncTools::syncPointList
        (
            mesh(),
            relaxationLevel,
            maxEqOp<label>(),
            label(0)
        );
    }

    // Check for convergence
    const label count(countPos(relaxationLevel));
    const bool converged(count == 0);

    //if (converged)
    //{
    //    Info<< "... Converged" << endl << endl;
    //}
    //else
    //{
    //    Info<< "... Not converged" << endl << endl;
    //}

    return converged;
}
Foam::label Foam::hexMeshSmootherMotionSolver::countPos
(
    const labelList& elems
) const
{
    label n = 0;
    for (const label elem : elems)
    {
        if (elem > 0)
        {
            n++;
        }
    }
    return returnReduce(n, sumOp<label>());
}


Foam::labelList Foam::hexMeshSmootherMotionSolver::countZeroOrPos
(
    const label size,
    const labelList& elems
) const
{
    labelList n(size, 0);
    for (const label elem : elems)
    {
        if (elem >= 0)
        {
            n[elem]++;
        }
    }
    Pstream::listCombineGather(n, plusEqOp<label>());
    Pstream::broadcast(n);
    return n;
}


void Foam::hexMeshSmootherMotionSolver::select
(
    const labelUList& lst,
    const label val,
    bitSet& isVal
) const
{
    isVal.set(lst.size());
    isVal = false;
    forAll(lst, i)
    {
        isVal[i] = (lst[i] == val);
    }
}


void Foam::hexMeshSmootherMotionSolver::laplaceSmooth
(
    const label type,
    const pointField& initialPoints,
    pointField& newPoints
) const
{
    if (initialPoints.size() != mesh().nPoints())
    {
        FatalErrorInFunction << "mesh().nPoints:" << mesh().nPoints()
            << " initial:" << initialPoints.size() << exit(FatalError);
    }

    newPoints.setSize(initialPoints.size());
    newPoints = Zero;
    labelList n(initialPoints.size(), 0);

    DynamicList<label> storage;
    forAll(pointTypes_, pointi)
    {
        if (pointTypes_[pointi] == INTERIOR)
        {
            const labelList& pPoints = mesh().pointPoints(pointi, storage);
            for (const label otherPointi : pPoints)
            {
                if (isMasterPoint_[otherPointi])
                {
                    newPoints[pointi] += initialPoints[otherPointi];
                    n[pointi]++;
                }
            }
            //Pout<< "Moving internal point " << initialPoints[pointi]
            //    << " to average " << newPoints[pointi]/n[pointi]
            //    << " of " << n[pointi] << " points" << endl;
        }
    }

    // Combine
    syncTools::syncPointList
    (
        mesh(),
        n,
        plusEqOp<label>(),
        label(0)
    );
    syncTools::syncPointList
    (
        mesh(),
        newPoints,
        plusEqOp<vector>(),
        vector::zero
    );
    forAll(newPoints, pointi)
    {
        if (n[pointi] == 0)
        {
            // This can happen if not interior point
            newPoints[pointi] = initialPoints[pointi];
            //Pout<< "Not Moving boundary point " << newPoints[pointi] << endl;
        }
        else
        {
            newPoints[pointi] /= n[pointi];
            //Pout<< "Moving internal point " << initialPoints[pointi]
            //    << " to " << newPoints[pointi] << endl;
        }
    }
}
void Foam::hexMeshSmootherMotionSolver::featLaplaceSmooth
(
    const indirectPrimitivePatch& pp,
    const pointField& initialPoints,
    pointField& newPoints
) const
{
    if (initialPoints.size() != pp.nPoints())
    {
        FatalErrorInFunction << "pp.nPoints:" << pp.nPoints()
            << " initial:" << initialPoints.size() << exit(FatalError);
    }

    newPoints.setSize(pp.nPoints());
    newPoints = Zero;
    labelList n(pp.nPoints(), 0);

    const edgeList& edges = pp.edges();
    const labelListList& pointEdges = pp.pointEdges();
    const labelList& meshPoints = pp.meshPoints();

    forAll(pointEdges, pointi)
    {
        const label myConstraint = pointTypes_[meshPoints[pointi]];
        if (myConstraint != INTERIOR)   // pp points should never be interior
        {
            const labelList& pEdges = pointEdges[pointi];
            //Pout<< "For boundary point:" << initialPoints[pointi]
            //    << endl;

            for (const label edgei : pEdges)
            {
                const label otherPointi = edges[edgei].otherVertex(pointi);
                const label otherMeshPointi = meshPoints[otherPointi];
                const label otherConstraint = pointTypes_[otherMeshPointi];

                if
                (
                    (otherConstraint != INTERIOR)   // Should not happen
                 && (myConstraint <= otherConstraint)
                 && isMasterPoint_[otherMeshPointi]
                )
                {
                    //Pout<< "    summing boundary point:"
                    //    << initialPoints[otherPointi] << endl;

                    newPoints[pointi] += initialPoints[otherPointi];
                    n[pointi]++;
                }
            }
        }
    }

    // Combine
    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        n,
        plusEqOp<label>(),
        label(0)
    );
    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        newPoints,
        plusEqOp<vector>(),
        vector::zero
    );

    forAll(newPoints, pointi)
    {
        if (n[pointi] == 0)
        {
            // This can happen if surface point surrounded by feature points
            // only.
            newPoints[pointi] = initialPoints[pointi];
            //Pout<< "Not Moving boundary point " << newPoints[pointi] << endl;
        }
        else
        {
            newPoints[pointi] /= n[pointi];
            //Pout<< "Moving surface point " << initialPoints[pointi]
            //    << " to average " << newPoints[pointi]
            //    << " of " << n[pointi] << " points" << endl;
        }
    }
}
void Foam::hexMeshSmootherMotionSolver::snapBoundaryPoints
(
    const scalar scale,
    const pointField& initialPoints,
    pointField& newPoints
) const
{
    if (initialPoints.size() != pointDisplacement_.mesh().size())
    {
        FatalErrorInFunction
            << "mesh.nPoints():" << pointDisplacement_.mesh().size()
            << " initial:" << initialPoints.size() << exit(FatalError);
    }
    const indirectPrimitivePatch& bnd0 = bnd0Ptr_();
    const labelList& mp = bnd0.meshPoints();

    // Save old point location
    const vectorField bndPoints(initialPoints, mp);

    // Update pointDisplacement_ to be consistent with mesh points being set to
    // initialPoints. This makes sure that the snapping is done using the
    // initialPoints as starting point
    pointDisplacement_.primitiveFieldRef() = initialPoints-points0();
    // 'snap' using boundary conditions
    pointDisplacement_.correctBoundaryConditions();

    // Calculate new position
    newPoints = points0() + pointDisplacement().internalField();

    if (scale < 1.0)
    {
        // Underrelax
        vectorField d(newPoints, mp);
        d -= bndPoints;
        d *= scale;
        d += bndPoints;
        UIndirectList<point>(newPoints, mp) = d;
    }
}

//void Foam::hexMeshSmootherMotionSolver::writeOBJ
//(
//    const fileName& name,
//    const pointField& p0,
//    const pointField& p1
//) const
//{
//    OBJstream os(mesh().time().path()/name);
//    forAll(p0, pointi)
//    {
//        os.write(linePointRef(p0[pointi], p1[pointi]));
//    }
//    Pout<< "Dumped to " << os.name() << endl;
//}
//void Foam::hexMeshSmootherMotionSolver::writeOBJ
//(
//    const fileName& name,
//    const UIndirectList<face>& pp,
//    const pointField& points
//) const
//{
//    const faceList fcs(pp);
//
//    OBJstream os(mesh().time().path()/name);
//    os.write(fcs, points, false);
//    Pout<< "Dumped faces to " << os.name() << endl;
//}


void Foam::hexMeshSmootherMotionSolver::emptyCorrectPoints
(
    pointVectorField& pointDisplacement
) const
{
    // Assume empty point patches are already in correct location
    // so knock out any off-plane displacement.
    auto& fld = pointDisplacement.primitiveFieldRef();
    for (const auto& ppf : pointDisplacement.boundaryField())
    {
        if (isA<emptyPointPatchVectorField>(ppf))
        {
            const auto& mp = ppf.patch().meshPoints();
            forAll(mp, i)
            {
                pointConstraint pc;
                ppf.patch().applyConstraint(i, pc);
                fld[mp[i]] = pc.constrainDisplacement(fld[mp[i]]);
            }
        }
    }

    pointField wantedPoints(points0() + fld);
    twoDCorrectPoints(wantedPoints);
    fld = wantedPoints-points0();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexMeshSmootherMotionSolver::
hexMeshSmootherMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    pointSmoother_(pointSmoother::New(mesh, coeffDict())),
    nPointSmootherIter_
    (
        readLabel(coeffDict().lookup("nPointSmootherIter"))
    ),
    relaxationFactors_(coeffDict().lookup("relaxationFactors")),
    relaxationLevel_(mesh.nPoints(), 0),
    relaxedPoints_(mesh.points()),
    //surfacesDict_(coeffDict().subDict("geometry")),
    //featureAngle_(coeffDict().get<scalar>("featureAngle")),
    //snapPatches_
    //(
    //    mesh.boundaryMesh().patchSet
    //    (
    //        coeffDict().get<wordRes>("patches")
    //    ).sortedToc()
    //),
    //snapZones_
    //(
    //    mesh.faceZones().indices
    //    (
    //        coeffDict().get<wordRes>("faceZones")
    //    )
    //),
    snapScale_(Function1<scalar>::New("snapScale", coeffDict())),
    isMasterPoint_(syncTools::getMasterPoints(mesh)),
    // Create big primitivePatch for all outside and any features on it
    //bnd0Ptr_(makePatch(mesh, snapPatches_, snapZones_, points0()))
    bnd0Ptr_
    (
        makePatch
        (
            mesh,
            nonConstraintPatches(mesh),
            labelList::null(),
            points0()
        )
    )
{
//     findSurfaces();
//
//     // Do multi-patch constraints
//     // ~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//     const scalar featureEdgeCos(Foam::cos(featureAngle_));
//     const scalar featurePointCos(featureEdgeCos);
//
//     const indirectPrimitivePatch& bnd0 = bnd0Ptr_();
//
//     calcConstraints
//     (
//         featureEdgeCos,
//         featurePointCos,
//         bnd0,
//         bnd0EdgeConstraints_,
//         bnd0PointConstraints_
//     );

    const pointMesh& pMesh = pointMesh::New(mesh, IOobject::READ_IF_PRESENT);

    pointTypes_.setSize(pMesh.size());
    pointTypes_ = INTERIOR;

    for (const auto& pp : pMesh.boundary())
    {
        if (!isA<meshPointPatch>(pp))
        {
            const auto& mp = pp.meshPoints();
            UIndirectList<label>(pointTypes_, mp) = pointType::SURFACE;
        }
    }

    // Override with any explicit constraint boundaries
    for (const auto& pp : pMesh.boundary())
    {
        const auto* meshPointPtr = isA<meshPointPatch>(pp);
        if (meshPointPtr)
        {
            const auto& constraints = meshPointPtr->constraints();
            const auto& mp = meshPointPtr->meshPoints();

            forAll(mp, i)
            {
                pointTypes_[mp[i]] = pointType(constraints[i].first());
            }
        }
    }

    // Make sure coupled points agree. Max constraint wins.
    syncTools::syncPointList
    (
        mesh,
        pointTypes_,
        maxEqOp<label>(),
        label(0)
    );

    bitSet isVal;
    select(pointTypes_, POINT, isVal);
    const label nFeatPoint = returnReduce(isVal.count(), sumOp<label>());
    select(pointTypes_, EDGE, isVal);
    const label nFeatEdge = returnReduce(isVal.count(), sumOp<label>());
    select(pointTypes_, SURFACE, isVal);
    const label nSurface = returnReduce(isVal.count(), sumOp<label>());
    select(pointTypes_, INTERIOR, isVal);
    const label nInternal = returnReduce(isVal.count(), sumOp<label>());
    Info<< "Attraction:" << nl
        << "    feature point:" << nFeatPoint << nl
        << "    feature edge :" << nFeatEdge << nl
        << "    surface      :" << nSurface << nl
        << "    none         :" << nInternal
        << endl;
}


Foam::hexMeshSmootherMotionSolver::
hexMeshSmootherMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName),
    pointSmoother_(pointSmoother::New(mesh, coeffDict())),
    //pointSmoother_
    //(
    //    pointSmoother::New
    //    (
    //        coeffDict(),
    //        displacementMotionSolver::pointDisplacement()
    //    )
    //),
    nPointSmootherIter_
    (
        readLabel(coeffDict().lookup("nPointSmootherIter"))
    ),
    relaxationFactors_(coeffDict().lookup("relaxationFactors")),
    relaxationLevel_(mesh.nPoints(), 0),
    relaxedPoints_(mesh.points()),
    //surfacesDict_(coeffDict().subDict("geometry")),
    //featureAngle_(coeffDict().get<scalar>("featureAngle")),
    //snapPatches_
    //(
    //    mesh.boundaryMesh().patchSet
    //    (
    //        coeffDict().get<wordReList>("patches")
    //    ).sortedToc()
    //),
    //snapZones_
    //(
    //    mesh.faceZones().indices
    //    (
    //        coeffDict().get<wordRes>("faceZones")
    //    )
    //),
    snapScale_(Function1<scalar>::New("snapScale", coeffDict())),
    isMasterPoint_(syncTools::getMasterPoints(mesh)),
    // Create big primitivePatch for all outside and any features on it
    //bnd0Ptr_(makePatch(mesh, snapPatches_, snapZones_, points0))
    bnd0Ptr_
    (
        makePatch
        (
            mesh,
            nonConstraintPatches(mesh),
            labelList::null(),
            points0
        )
    )
{
//    findSurfaces();
//
//    const scalar featureEdgeCos(Foam::cos(featureAngle_));
//    const scalar featurePointCos(featureEdgeCos);
//
//    const indirectPrimitivePatch& bnd0 = bnd0Ptr_();
//
//    calcConstraints
//    (
//        featureEdgeCos,
//        featurePointCos,
//        bnd0,
//        bnd0EdgeConstraints_,
//        bnd0PointConstraints_
//    );

    const pointMesh& pMesh = pointMesh::New(mesh, IOobject::READ_IF_PRESENT);

    pointTypes_.setSize(mesh.nPoints());
    pointTypes_ = INTERIOR;

    for (const auto& pp : pMesh.boundary())
    {
        if (!isA<meshPointPatch>(pp) && !pp.coupled())
        {
            const auto& mp = pp.meshPoints();
            UIndirectList<label>(pointTypes_, mp) = pointType::SURFACE;
        }
    }

    // Override with any explicit constraint boundaries
    for (const auto& pp : pMesh.boundary())
    {
        const auto* meshPointPtr = isA<meshPointPatch>(pp);
        if (meshPointPtr)
        {
            const auto& constraints = meshPointPtr->constraints();
            const auto& mp = meshPointPtr->meshPoints();

            forAll(mp, i)
            {
                pointTypes_[mp[i]] = pointType(constraints[i].first());
            }
        }
    }

    // Make sure coupled points agree. Max constraint wins.
    syncTools::syncPointList
    (
        mesh,
        pointTypes_,
        maxEqOp<label>(),
        label(0)
    );

    bitSet isVal;
    select(pointTypes_, POINT, isVal);
    const label nFeatPoint = returnReduce(isVal.count(), sumOp<label>());
    select(pointTypes_, EDGE, isVal);
    const label nFeatEdge = returnReduce(isVal.count(), sumOp<label>());
    select(pointTypes_, SURFACE, isVal);
    const label nSurface = returnReduce(isVal.count(), sumOp<label>());
    select(pointTypes_, INTERIOR, isVal);
    const label nInternal = returnReduce(isVal.count(), sumOp<label>());
    Info<< "Attraction:" << nl
        << "    feature point:" << nFeatPoint << nl
        << "    feature edge :" << nFeatEdge << nl
        << "    surface      :" << nSurface << nl
        << "    none         :" << nInternal
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hexMeshSmootherMotionSolver::
~hexMeshSmootherMotionSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::hexMeshSmootherMotionSolver::curPoints() const
{
    //Note: twoDCorrect already done by ::solve
    return relaxedPoints_;
}


void Foam::hexMeshSmootherMotionSolver::solve()
{
    // Update any internal storage for current state
    movePoints(mesh().points());

    // No updating of bc since we don't want to do snapping yet - is done
    // later on
    //pointDisplacement().boundaryFieldRef().updateCoeffs();

    const indirectPrimitivePatch& bnd0 = bnd0Ptr_();
    const labelList& mp = bnd0.meshPoints();

    // Points on boundary
    bitSet isBndPoint(mesh().nPoints(), false);
    isBndPoint.set(mp);

    // Points not on boundary
    bitSet isInternalPoint;
    select(pointTypes_, INTERIOR, isInternalPoint);


    const pointField& initialPoints = mesh().points();

    // Wanted locations - starts off from current mesh points. Note : could
    // use relaxedPoints_ but this should be equal to current mesh points unless
    // we have multiple motion solvers ...
    pointField movedPoints(initialPoints);


    // -1 indicates a point which is not to be moved
    //  0 is the starting value for a moving point
    const label nRelaxed = countPos(relaxationLevel_);
    //Pout<< "Starting relaxed:" << nRelaxed << endl;

    if (nRelaxed > 0)
    {
        // MeshSmoother::snapSmoothing()
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // - set initialPoint from relaxedPoint
        // - smooth internal points
        //      - and relax point motion (= adapt relaxationLevel_)
        // - snap boundary points
        //      - and relax point motion (= adapt relaxationLevel_)
        // - mark all points that are not snapped (i.e. snapped but relaxed)
        // - smooth all boundary points (using features)
        //      - and relax point motion (= adapt relaxationLevel_)

        // Laplace smoothing of interior (=unconstrained) points
        laplaceSmooth(INTERIOR, initialPoints, movedPoints);

        // Apply constraints
        //constrainDisplacement(movedPoints);

        // Determine relaxation level to move points
        relax
        (
            relaxationFactors_,
            isInternalPoint,
            initialPoints,      // starting location
            movedPoints,        // wanted location
            relaxedPoints_,     // inbetween location without errors
            relaxationLevel_
        );
        //Pout<< "After laplaceSmooth:"
        //    << countZeroOrPos(relaxationFactors_.size(), relaxationLevel_)
        //    << endl;
        //writeOBJ
        //(
        //    "laplaceSmooth_relax_" + mesh().time().timeName() + ".obj",
        //    bnd0,
        //    relaxedPoints_
        //);


        // Snap boundary points
        const scalar scale = snapScale_->value(mesh().time().timeIndex());
        if (scale > 0)
        {
            // snap to surface (=apply boundary conditions)
            snapBoundaryPoints(scale, initialPoints, movedPoints);

            // Apply constraints
            //constrainDisplacement(movedPoints);

            relax
            (
                relaxationFactors_,
                isBndPoint,
                initialPoints,
                movedPoints,
                relaxedPoints_,
                relaxationLevel_
            );
            //Pout<< "After snapping:"
            //    << countZeroOrPos(relaxationFactors_.size(), relaxationLevel_)
            //    << endl;
            //writeOBJ("snap_relax.obj", initialPoints, relaxedPoints_);

            // Now relaxedPoints_ with relaxationLevel_ 0 are perfectly snapped
            // (only applicable for bnd0 points)
        }

        // Laplace smoothing of (now snapped&relaxed) boundary points. Use
        // average of surrounding boundary points of same type only
        pointField bndMovedPoints;
        featLaplaceSmooth
        (
            bnd0,
            pointField(relaxedPoints_, mp),
            bndMovedPoints
        );
        UIndirectList<point>(movedPoints, mp) = bndMovedPoints;

        // Apply constraints
        //constrainDisplacement(movedPoints);
        //writeOBJ("featLaplaceSmooth.obj", initialPoints, movedPoints);

        relax
        (
            relaxationFactors_,
            isBndPoint,
            initialPoints,
            movedPoints,
            relaxedPoints_,
            relaxationLevel_
        );
    }
    else
    {
        // MeshSmoother::GETMeSmoothing()
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // - set initialPoint from relaxedPoint
        // - smooth all boundary points (using features)
        //      - and relax point motion (= adapt relaxationLevel_)
        // - snap boundary points
        //      - and relax point motion (= adapt relaxationLevel_)
        // - do all points GETMe
        //      - and relax point motion (= adapt relaxationLevel_)


        // Laplace smoothing of boundary points. Use average of
        // surrounding boundary points of same type only. TBD: make consistent
        // wrt coupled points not on boundary patch.
        pointField bndMovedPoints;
        featLaplaceSmooth
        (
            bnd0,
            pointField(relaxedPoints_, mp),
            bndMovedPoints
        );
        UIndirectList<point>(movedPoints, mp) = bndMovedPoints;
        //writeOBJ
        //(
        //    "featLaplaceSmooth_unconstrain_"
        //  + mesh().time().timeName()
        //  + ".obj",
        //    pointField(initialPoints, mp),
        //    bndMovedPoints
        //);

        // Apply constraints
        //constrainDisplacement(movedPoints);

        relax
        (
            relaxationFactors_,
            isBndPoint,
            initialPoints,
            movedPoints,
            relaxedPoints_,
            relaxationLevel_
        );
        //writeOBJ
        //(
        //    "featLaplaceSmooth_relax"
        //  + mesh().time().timeName()
        //  + ".obj",
        //    initialPoints,
        //    relaxedPoints_
        //);

        // Snap boundary points
        const scalar scale = snapScale_->value(mesh().time().timeIndex());
        if (scale > 0)
        {
            // snap to surface (=apply boundary conditions)
            snapBoundaryPoints(scale, relaxedPoints_, movedPoints);

            // Apply constraints
            //constrainDisplacement(movedPoints);
            //writeOBJ("snap.obj", initialPoints, movedPoints);

            relax
            (
                relaxationFactors_,
                isBndPoint,
                initialPoints,
                movedPoints,
                relaxedPoints_,
                relaxationLevel_
            );
            //writeOBJ("snap_relax.obj", initialPoints, relaxedPoints_);
        }

        vectorField fCtrs(mesh().nFaces());
        vectorField fAreas(mesh().nFaces());
        vectorField cellCtrs(mesh().nCells());
        scalarField cellVols(mesh().nCells());

        movedPoints = relaxedPoints_;

        for(label i = 0; i < nPointSmootherIter_; i ++)
        {
            // Starting from current points do smoothing
            // Calculate mesh quantities with new locations

            const auto& geom =
                reinterpret_cast<const fvMesh&>(mesh()).geometry();
            geom.updateGeom
            (
                movedPoints,
                mesh().points(),    // old points
                fCtrs,
                fAreas,
                cellCtrs,
                cellVols
            );

            //- Smooth point positions (returned as pointDisplacement w.r.t.
            //  points0)
            pointSmoother_->update
            (
                identity(mesh().nFaces()),
                points0(),
                movedPoints,
                fCtrs,
                fAreas,
                cellCtrs,
                cellVols,
                pointDisplacement_
            );
            // Keep points on empty patches
            emptyCorrectPoints(pointDisplacement());

            // Update moving points
            movedPoints = points0() + pointDisplacement().internalField();

            //// snap to surface (=apply boundary conditions)
            //snapBoundaryPoints(scale, movedPoints, movedPoints);
        }

        // snap to surface (=apply boundary conditions)
        //snapBoundaryPoints(scale, movedPoints, movedPoints);

        //writeOBJ
        //(
        //    "GETMeSmoothing_snapped_"
        //  + mesh().time().timeName()
        //  + ".obj",
        //    pointField(initialPoints, mp),
        //    pointField(movedPoints, mp)
        //);

        relax
        (
            relaxationFactors_,
            bitSet(mesh().nPoints(), true),
            initialPoints,
            movedPoints,
            relaxedPoints_,
            relaxationLevel_
        );
        //writeOBJ("GETMeSmoothing_relax.obj", initialPoints, relaxedPoints_);
    }
}


// ************************************************************************* //

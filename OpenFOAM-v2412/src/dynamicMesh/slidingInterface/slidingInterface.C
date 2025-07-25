/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "slidingInterface.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "plane.H"

// Index of debug signs:
// p - adjusting a projection point
// * - adjusting edge intersection

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(slidingInterface, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        slidingInterface,
        dictionary
    );
}


const Foam::Enum
<
    Foam::slidingInterface::typeOfMatch
>
Foam::slidingInterface::typeOfMatchNames
({
    { typeOfMatch::INTEGRAL, "integral" },
    { typeOfMatch::PARTIAL, "partial" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::checkDefinition()
{
    const polyMesh& mesh = topoChanger().mesh();

    if
    (
        !masterFaceZoneID_.active()
     || !slaveFaceZoneID_.active()
     || !cutPointZoneID_.active()
     || !cutFaceZoneID_.active()
     || !masterPatchID_.active()
     || !slavePatchID_.active()
    )
    {
        FatalErrorInFunction
            << "Not all zones and patches needed in the definition "
            << "have been found.  Please check your mesh definition."
            << abort(FatalError);
    }

    // Check the sizes and set up state
    if
    (
        mesh.faceZones()[masterFaceZoneID_.index()].empty()
     || mesh.faceZones()[slaveFaceZoneID_.index()].empty()
    )
    {
        FatalErrorInFunction
            << "Please check your mesh definition."
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "Sliding interface object " << name() << " :" << nl
            << "    master face zone: " << masterFaceZoneID_.index() << nl
            << "    slave face zone: " << slaveFaceZoneID_.index() << endl;
    }
}


void Foam::slidingInterface::clearOut() const
{
    clearPointProjection();
    clearAttachedAddressing();
    clearAddressing();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slidingInterface::slidingInterface
(
    const word& name,
    const label index,
    const polyTopoChanger& mme,
    const word& masterFaceZoneName,
    const word& slaveFaceZoneName,
    const word& cutPointZoneName,
    const word& cutFaceZoneName,
    const word& masterPatchName,
    const word& slavePatchName,
    const typeOfMatch tom,
    const bool coupleDecouple,
    const intersection::algorithm algo
)
:
    polyMeshModifier(name, index, mme, true),
    masterFaceZoneID_
    (
        masterFaceZoneName,
        mme.mesh().faceZones()
    ),
    slaveFaceZoneID_
    (
        slaveFaceZoneName,
        mme.mesh().faceZones()
    ),
    cutPointZoneID_
    (
        cutPointZoneName,
        mme.mesh().pointZones()
    ),
    cutFaceZoneID_
    (
        cutFaceZoneName,
        mme.mesh().faceZones()
    ),
    masterPatchID_
    (
        masterPatchName,
        mme.mesh().boundaryMesh()
    ),
    slavePatchID_
    (
        slavePatchName,
        mme.mesh().boundaryMesh()
    ),
    matchType_(tom),
    coupleDecouple_(coupleDecouple),
    attached_(false),
    projectionAlgo_(algo),
    trigger_(false),
    pointMergeTol_(pointMergeTolDefault_),
    edgeMergeTol_(edgeMergeTolDefault_),
    nFacesPerSlaveEdge_(nFacesPerSlaveEdgeDefault_),
    edgeFaceEscapeLimit_(edgeFaceEscapeLimitDefault_),
    integralAdjTol_(integralAdjTolDefault_),
    edgeMasterCatchFraction_(edgeMasterCatchFractionDefault_),
    edgeCoPlanarTol_(edgeCoPlanarTolDefault_),
    edgeEndCutoffTol_(edgeEndCutoffTolDefault_),
    cutFaceMasterPtr_(nullptr),
    cutFaceSlavePtr_(nullptr),
    masterFaceCellsPtr_(nullptr),
    slaveFaceCellsPtr_(nullptr),
    masterStickOutFacesPtr_(nullptr),
    slaveStickOutFacesPtr_(nullptr),
    retiredPointMapPtr_(nullptr),
    cutPointEdgePairMapPtr_(nullptr),
    slavePointPointHitsPtr_(nullptr),
    slavePointEdgeHitsPtr_(nullptr),
    slavePointFaceHitsPtr_(nullptr),
    masterPointEdgeHitsPtr_(nullptr),
    projectedSlavePointsPtr_(nullptr)
{
    checkDefinition();

    if (attached_)
    {
        FatalErrorInFunction
            << "Creation of a sliding interface from components "
            << "in attached state not supported."
            << abort(FatalError);
    }
    else
    {
        calcAttachedAddressing();
    }
}


Foam::slidingInterface::slidingInterface
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, dict.get<bool>("active")),
    masterFaceZoneID_
    (
        dict.get<keyType>("masterFaceZoneName"),
        mme.mesh().faceZones()
    ),
    slaveFaceZoneID_
    (
        dict.get<keyType>("slaveFaceZoneName"),
        mme.mesh().faceZones()
    ),
    cutPointZoneID_
    (
        dict.get<keyType>("cutPointZoneName"),
        mme.mesh().pointZones()
    ),
    cutFaceZoneID_
    (
        dict.get<keyType>("cutFaceZoneName"),
        mme.mesh().faceZones()
    ),
    masterPatchID_
    (
        dict.get<keyType>("masterPatchName"),
        mme.mesh().boundaryMesh()
    ),
    slavePatchID_
    (
        dict.get<keyType>("slavePatchName"),
        mme.mesh().boundaryMesh()
    ),
    matchType_(typeOfMatchNames.get("typeOfMatch", dict)),
    coupleDecouple_(dict.get<bool>("coupleDecouple")),
    attached_(dict.get<bool>("attached")),
    projectionAlgo_
    (
        intersection::algorithmNames_.get("projection", dict)
    ),
    trigger_(false),
    cutFaceMasterPtr_(nullptr),
    cutFaceSlavePtr_(nullptr),
    masterFaceCellsPtr_(nullptr),
    slaveFaceCellsPtr_(nullptr),
    masterStickOutFacesPtr_(nullptr),
    slaveStickOutFacesPtr_(nullptr),
    retiredPointMapPtr_(nullptr),
    cutPointEdgePairMapPtr_(nullptr),
    slavePointPointHitsPtr_(nullptr),
    slavePointEdgeHitsPtr_(nullptr),
    slavePointFaceHitsPtr_(nullptr),
    masterPointEdgeHitsPtr_(nullptr),
    projectedSlavePointsPtr_(nullptr)
{
    // Optionally default tolerances from dictionary
    setTolerances(dict);

    checkDefinition();

    // If the interface is attached, the master and slave face zone addressing
    // needs to be read in; otherwise it will be created
    if (attached_)
    {
        if (debug)
        {
            Pout<< "slidingInterface::slidingInterface(...) "
                << " for object " << name << " : "
                << "Interface attached.  Reading master and slave face zones "
                << "and retired point lookup." << endl;
        }

        // The face zone addressing is written out in the definition dictionary
        masterFaceCellsPtr_.reset(new labelList());
        slaveFaceCellsPtr_.reset(new labelList());
        masterStickOutFacesPtr_.reset(new labelList());
        slaveStickOutFacesPtr_.reset(new labelList());
        retiredPointMapPtr_.reset(new Map<label>());
        cutPointEdgePairMapPtr_.reset(new Map<Pair<edge>>());

        dict.readEntry("masterFaceCells", *masterFaceCellsPtr_);
        dict.readEntry("slaveFaceCells", *slaveFaceCellsPtr_);
        dict.readEntry("masterStickOutFaces", *masterStickOutFacesPtr_);
        dict.readEntry("slaveStickOutFaces", *slaveStickOutFacesPtr_);
        dict.readEntry("retiredPointMap", *retiredPointMapPtr_);
        dict.readEntry("cutPointEdgePairMap", *cutPointEdgePairMapPtr_);
    }
    else
    {
        calcAttachedAddressing();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slidingInterface::clearAddressing() const
{
    cutFaceMasterPtr_.reset(nullptr);
    cutFaceSlavePtr_.reset(nullptr);
}


const Foam::faceZoneID& Foam::slidingInterface::masterFaceZoneID() const
{
    return masterFaceZoneID_;
}


const Foam::faceZoneID& Foam::slidingInterface::slaveFaceZoneID() const
{
    return slaveFaceZoneID_;
}


bool Foam::slidingInterface::changeTopology() const
{
    if (coupleDecouple_)
    {
        // Always changes.  If not attached, project points
        if (debug)
        {
            Pout<< "bool slidingInterface::changeTopology() const "
                << "for object " << name() << " : "
                << "Couple-decouple mode." << endl;
        }

        if (!attached_)
        {
            projectPoints();
        }
        else
        {
        }

        return true;
    }

    if
    (
        attached_
     && !topoChanger().mesh().changing()
    )
    {
        // If the mesh is not moving or morphing and the interface is
        // already attached, the topology will not change
        return false;
    }
    else
    {
        // Check if the motion changes point projection
        return projectPoints();
    }
}


void Foam::slidingInterface::setRefinement(polyTopoChange& ref) const
{
    if (coupleDecouple_)
    {
        if (attached_)
        {
            // Attached, coupling
            decoupleInterface(ref);
        }
        else
        {
            // Detached, coupling
            coupleInterface(ref);
        }

        return;
    }

    if (trigger_)
    {
        if (attached_)
        {
            // Clear old coupling data
            clearCouple(ref);
        }

        coupleInterface(ref);

        trigger_ = false;
    }
}


void Foam::slidingInterface::modifyMotionPoints(pointField& motionPoints) const
{
    if (debug)
    {
        Pout<< "void slidingInterface::modifyMotionPoints("
            << "pointField& motionPoints) const for object " << name() << " : "
            << "Adjusting motion points." << endl;
    }

    const polyMesh& mesh = topoChanger().mesh();

    // Get point from the cut zone
    const labelList& cutPoints = mesh.pointZones()[cutPointZoneID_.index()];

    if (cutPoints.size() && !projectedSlavePointsPtr_)
    {
        return;
    }
    else
    {
        const pointField& projectedSlavePoints = *projectedSlavePointsPtr_;

        const Map<label>& rpm = retiredPointMap();

        const Map<Pair<edge>>& cpepm = cutPointEdgePairMap();

        const Map<label>& slaveZonePointMap =
            mesh.faceZones()[slaveFaceZoneID_.index()]().meshPointMap();

        const primitiveFacePatch& masterPatch =
            mesh.faceZones()[masterFaceZoneID_.index()]();
        const edgeList& masterEdges = masterPatch.edges();
        const pointField& masterLocalPoints = masterPatch.localPoints();

        const primitiveFacePatch& slavePatch =
            mesh.faceZones()[slaveFaceZoneID_.index()]();
        const edgeList& slaveEdges = slavePatch.edges();
        const pointField& slaveLocalPoints = slavePatch.localPoints();
        const vectorField& slavePointNormals = slavePatch.pointNormals();

        for (const label pointi : cutPoints)
        {
            // Try to find the cut point in retired points
            const auto rpmIter = rpm.cfind(pointi);

            if (rpmIter.good())
            {
                if (debug)
                {
                    Pout<< "p";
                }

                // Cut point is a retired point
                motionPoints[cutPoints[pointi]] =
                    projectedSlavePoints[slaveZonePointMap.find(rpmIter())()];
            }
            else
            {
                // A cut point is not a projected slave point.  Therefore, it
                // must be an edge-to-edge intersection.

                const auto cpepmIter = cpepm.cfind(pointi);

                if (cpepmIter.good())
                {
                    // Pout<< "Need to re-create hit for point "
                    //     << cutPoints[pointi]
                    //     << " lookup: " << cpepmIter()
                    //     << endl;

                    // Note.
                    // The edge cutting code is repeated in
                    // slidingInterface::coupleInterface.  This is done for
                    // efficiency reasons and avoids multiple creation of
                    // cutting planes.  Please update both simultaneously.
                    //
                    const edge& globalMasterEdge = cpepmIter().first();

                    const label curMasterEdgeIndex =
                        masterPatch.whichEdge
                        (
                            edge
                            (
                                masterPatch.whichPoint
                                (
                                    globalMasterEdge.start()
                                ),
                                masterPatch.whichPoint
                                (
                                    globalMasterEdge.end()
                                )
                            )
                        );

                    const edge& cme = masterEdges[curMasterEdgeIndex];

                    // Pout<< "curMasterEdgeIndex: " << curMasterEdgeIndex
                    //     << " cme: " << cme
                    //     << endl;

                    const edge& globalSlaveEdge = cpepmIter().second();

                    const label curSlaveEdgeIndex =
                        slavePatch.whichEdge
                        (
                            edge
                            (
                                slavePatch.whichPoint
                                (
                                    globalSlaveEdge.start()
                                ),
                                slavePatch.whichPoint
                                (
                                    globalSlaveEdge.end()
                                )
                            )
                        );

                    const edge& curSlaveEdge = slaveEdges[curSlaveEdgeIndex];
                    // Pout<< "curSlaveEdgeIndex: " << curSlaveEdgeIndex
                    //     << " curSlaveEdge: " << curSlaveEdge
                    //     << endl;
                    const point& a = projectedSlavePoints[curSlaveEdge.start()];
                    const point& b = projectedSlavePoints[curSlaveEdge.end()];

                    point c =
                        0.5*
                        (
                            slaveLocalPoints[curSlaveEdge.start()]
                          + slavePointNormals[curSlaveEdge.start()]
                          + slaveLocalPoints[curSlaveEdge.end()]
                          + slavePointNormals[curSlaveEdge.end()]
                        );

                    // Create the plane
                    plane cutPlane(a, b, c);

                    linePointRef curSlaveLine =
                        curSlaveEdge.line(slaveLocalPoints);
                    const scalar curSlaveLineMag = curSlaveLine.mag();

                    scalar cutOnMaster =
                        cutPlane.lineIntersect
                        (
                            cme.line(masterLocalPoints)
                        );

                    if
                    (
                        cutOnMaster > edgeEndCutoffTol_
                     && cutOnMaster < 1.0 - edgeEndCutoffTol_
                    )
                    {
                        // Master is cut, check the slave
                        point masterCutPoint =
                            masterLocalPoints[cme.start()]
                          + cutOnMaster*cme.vec(masterLocalPoints);

                        pointHit slaveCut =
                            curSlaveLine.nearestDist(masterCutPoint);

                        if (slaveCut.hit())
                        {
                            // Strict checking of slave cut to avoid capturing
                            // end points.
                            scalar cutOnSlave =
                                (
                                    (
                                        slaveCut.point()
                                      - curSlaveLine.start()
                                    ) & curSlaveLine.vec()
                                )/sqr(curSlaveLineMag);

                            // Calculate merge tolerance from the
                            // target edge length
                            scalar mergeTol =
                                edgeCoPlanarTol_*mag(b - a);

                            if
                            (
                                cutOnSlave > edgeEndCutoffTol_
                             && cutOnSlave < 1.0 - edgeEndCutoffTol_
                             && slaveCut.distance() < mergeTol
                            )
                            {
                                // Cut both master and slave.
                                motionPoints[pointi] = masterCutPoint;
                            }
                        }
                        else
                        {
                            Pout<< "Missed slave edge!!!  This is an error.  "
                                << "Master edge: "
                                << cme.line(masterLocalPoints)
                                << " slave edge: " << curSlaveLine
                                << " point: " << masterCutPoint
                                << " weight: " <<
                                (
                                    (
                                        slaveCut.point()
                                      - curSlaveLine.start()
                                    ) & curSlaveLine.vec()
                                )/sqr(curSlaveLineMag)
                                << endl;
                        }
                    }
                    else
                    {
                        Pout<< "Missed master edge!!!  This is an error"
                            << endl;
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << "Cut point " << cutPoints[pointi]
                        << " not recognised as either the projected "
                        << "or as intersection point.  Error in point "
                        << "projection or data mapping"
                        << abort(FatalError);
                }
            }
        }
        if (debug)
        {
            Pout<< endl;
        }
    }
}


void Foam::slidingInterface::updateMesh(const mapPolyMesh& m)
{
    if (debug)
    {
        Pout<< "void slidingInterface::updateMesh(const mapPolyMesh& m)"
            << " const for object " << name() << " : "
            << "Updating topology." << endl;
    }

    // Mesh has changed topologically.  Update local topological data
    const polyMesh& mesh = topoChanger().mesh();

    masterFaceZoneID_.update(mesh.faceZones());
    slaveFaceZoneID_.update(mesh.faceZones());
    cutPointZoneID_.update(mesh.pointZones());
    cutFaceZoneID_.update(mesh.faceZones());

    masterPatchID_.update(mesh.boundaryMesh());
    slavePatchID_.update(mesh.boundaryMesh());

//MJ:Disabled updating
//    if (!attached())
//    {
//        calcAttachedAddressing();
//    }
//    else
//    {
//        renumberAttachedAddressing(m);
//    }
}


const Foam::pointField& Foam::slidingInterface::pointProjection() const
{
    if (!projectedSlavePointsPtr_)
    {
        projectPoints();
    }

    return *projectedSlavePointsPtr_;
}


void Foam::slidingInterface::setTolerances(const dictionary&dict, bool report)
{
    pointMergeTol_ = dict.getOrDefault<scalar>
    (
        "pointMergeTol",
        pointMergeTol_
    );
    edgeMergeTol_ = dict.getOrDefault<scalar>
    (
        "edgeMergeTol",
        edgeMergeTol_
    );
    nFacesPerSlaveEdge_ = dict.getOrDefault<label>
    (
        "nFacesPerSlaveEdge",
        nFacesPerSlaveEdge_
    );
    edgeFaceEscapeLimit_ = dict.getOrDefault<label>
    (
        "edgeFaceEscapeLimit",
        edgeFaceEscapeLimit_
    );
    integralAdjTol_ = dict.getOrDefault<scalar>
    (
        "integralAdjTol",
        integralAdjTol_
    );
    edgeMasterCatchFraction_ = dict.getOrDefault<scalar>
    (
        "edgeMasterCatchFraction",
        edgeMasterCatchFraction_
    );
    edgeCoPlanarTol_ = dict.getOrDefault<scalar>
    (
        "edgeCoPlanarTol",
        edgeCoPlanarTol_
    );
    edgeEndCutoffTol_ = dict.getOrDefault<scalar>
    (
        "edgeEndCutoffTol",
        edgeEndCutoffTol_
    );

    if (report)
    {
        Info<< "Sliding interface parameters:" << nl
            << "pointMergeTol            : " << pointMergeTol_ << nl
            << "edgeMergeTol             : " << edgeMergeTol_ << nl
            << "nFacesPerSlaveEdge       : " << nFacesPerSlaveEdge_ << nl
            << "edgeFaceEscapeLimit      : " << edgeFaceEscapeLimit_ << nl
            << "integralAdjTol           : " << integralAdjTol_ << nl
            << "edgeMasterCatchFraction  : " << edgeMasterCatchFraction_ << nl
            << "edgeCoPlanarTol          : " << edgeCoPlanarTol_ << nl
            << "edgeEndCutoffTol         : " << edgeEndCutoffTol_ << endl;
    }
}


void Foam::slidingInterface::write(Ostream& os) const
{
    os  << nl << type() << nl
        << name()<< nl
        << masterFaceZoneID_.name() << nl
        << slaveFaceZoneID_.name() << nl
        << cutPointZoneID_.name() << nl
        << cutFaceZoneID_.name() << nl
        << masterPatchID_.name() << nl
        << slavePatchID_.name() << nl
        << typeOfMatchNames[matchType_] << nl
        << coupleDecouple_ << nl
        << attached_ << endl;
}


// To write out all those tolerances
#undef  WRITE_NON_DEFAULT
#define WRITE_NON_DEFAULT(name) \
    if ( name ## _ != name ## Default_ ) os.writeEntry( #name , name ## _ );


void Foam::slidingInterface::writeDict(Ostream& os) const
{
    os  << nl;

    os.beginBlock(name());

    os.writeEntry("type", type());
    os.writeEntry("masterFaceZoneName", masterFaceZoneID_.name());
    os.writeEntry("slaveFaceZoneName", slaveFaceZoneID_.name());
    os.writeEntry("cutPointZoneName", cutPointZoneID_.name());
    os.writeEntry("cutFaceZoneName", cutFaceZoneID_.name());
    os.writeEntry("masterPatchName", masterPatchID_.name());
    os.writeEntry("slavePatchName", slavePatchID_.name());
    os.writeEntry("typeOfMatch", typeOfMatchNames[matchType_]);
    os.writeEntry("coupleDecouple", coupleDecouple_);
    os.writeEntry("projection", intersection::algorithmNames_[projectionAlgo_]);
    os.writeEntry("attached", attached_);
    os.writeEntry("active", active());

    if (attached_)
    {
        masterFaceCellsPtr_->writeEntry("masterFaceCells", os);
        slaveFaceCellsPtr_->writeEntry("slaveFaceCells", os);
        masterStickOutFacesPtr_->writeEntry("masterStickOutFaces", os);
        slaveStickOutFacesPtr_->writeEntry("slaveStickOutFaces", os);

        os.writeEntry("retiredPointMap", retiredPointMap());
        os.writeEntry("cutPointEdgePairMap", cutPointEdgePairMap());
    }

    WRITE_NON_DEFAULT(pointMergeTol)
    WRITE_NON_DEFAULT(edgeMergeTol)
    WRITE_NON_DEFAULT(nFacesPerSlaveEdge)
    WRITE_NON_DEFAULT(edgeFaceEscapeLimit)
    WRITE_NON_DEFAULT(integralAdjTol)
    WRITE_NON_DEFAULT(edgeMasterCatchFraction)
    WRITE_NON_DEFAULT(edgeCoPlanarTol)
    WRITE_NON_DEFAULT(edgeEndCutoffTol)

    #undef WRITE_NON_DEFAULT

    os.endBlock();
}


// ************************************************************************* //

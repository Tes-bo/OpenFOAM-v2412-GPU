/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "wallBoundedParticle.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackCloudType>
void Foam::wallBoundedParticle::patchInteraction
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackFraction
)
{
    typename TrackCloudType::particleType& p =
        static_cast<typename TrackCloudType::particleType&>(*this);
    typename TrackCloudType::particleType::trackingData& ttd =
        static_cast<typename TrackCloudType::particleType::trackingData&>(td);

    if (!mesh().isInternalFace(face()))
    {
        label origFacei = face();
        label patchi = patch();

        // Did patch interaction model switch patches?
        // Note: recalculate meshEdgeStart_, diagEdge_!
        if (face() != origFacei)
        {
            patchi = patch();
        }

        const polyPatch& patch = mesh().boundaryMesh()[patchi];

        if (isA<processorPolyPatch>(patch))
        {
            p.hitProcessorPatch(cloud, ttd);
        }
        else if (isA<wallPolyPatch>(patch))
        {
            p.hitWallPatch(cloud, ttd);
        }
        else
        {
            td.keepParticle = false;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackCloudType>
Foam::scalar Foam::wallBoundedParticle::trackToEdge
(
    TrackCloudType& cloud,
    trackingData& td,
    const vector& endPosition
)
{
    // Track particle to a given position and returns 1.0 if the
    // trajectory is completed without hitting a face otherwise
    // stops at the face and returns the fraction of the trajectory
    // completed.
    // on entry 'stepFraction()' should be set to the fraction of the
    // time-step at which the tracking starts.

    // Are we on a track face? If not we do a topological walk.

    // Particle:
    // - cell_              always set
    // - tetFace_, tetPt_   always set (these identify tet particle is in)
    // - optionally meshEdgeStart_ or  diagEdge_ set (edge particle is on)

    //checkInside();
    //checkOnTriangle(localPosition_);
    //if (meshEdgeStart_ != -1 || diagEdge_ != -1)
    //{
    //    checkOnEdge();
    //}

    scalar trackFraction = 0.0;

    if (!td.isWallPatch_[tetFace()])
    {
        // Don't track across face. Just walk in cell. Particle is on
        // mesh edge (as indicated by meshEdgeStart_).
        const edge meshEdge(currentEdge());

        // If internal face check whether to go to neighbour cell or just
        // check to the other internal tet on the edge.
        if (mesh().isInternalFace(tetFace()))
        {
            label nbrCelli =
            (
                this->cell() == mesh().faceOwner()[face()]
              ? mesh().faceNeighbour()[face()]
              : mesh().faceOwner()[face()]
            );
            // Check angle to nbrCell tet. Is it in the direction of the
            // endposition? i.e. since volume of nbr tet is positive the
            // tracking direction should be into the tet.
            tetIndices nbrTi(nbrCelli, tetFace(), tetPt());

            const bool posVol = (nbrTi.tet(mesh()).mag() > 0);
            const vector path(endPosition - localPosition_);

            if (posVol == ((nbrTi.faceTri(mesh()).areaNormal() & path) < 0))
            {
                // Change into nbrCell. No need to change tetFace, tetPt.
                //Pout<< "    crossed from cell:" << celli_
                //    << " into " << nbrCelli << endl;
                this->cell() = nbrCelli;
                patchInteraction(cloud, td, trackFraction);
            }
            else
            {
                // Walk to other face on edge. Changes tetFace, tetPt but not
                // cell.
                crossEdgeConnectedFace(meshEdge);
                patchInteraction(cloud, td, trackFraction);
            }
        }
        else
        {
            // Walk to other face on edge. This might give loop since
            // particle should have been removed?
            crossEdgeConnectedFace(meshEdge);
            patchInteraction(cloud, td, trackFraction);
        }
    }
    else
    {
        // We're inside a tet on the wall. Check if the current tet is
        // the one to cross. If not we cross into the neighbouring triangle.

        if (mesh().isInternalFace(tetFace()))
        {
            FatalErrorInFunction
                << "Can only track on boundary faces."
                << " Face:" << tetFace()
                << " at:" << mesh().faceCentres()[tetFace()]
                << abort(FatalError);
        }

        const triFace tri(currentTetIndices().faceTriIs(mesh(), false));
        const vector n = tri.unitNormal(mesh().points());

        point projectedEndPosition = endPosition;

        const bool posVol = (currentTetIndices().tet(mesh()).mag() > 0);

        if (!posVol)
        {
            // Negative tet volume. Track back by setting the end point
            projectedEndPosition =
                localPosition_ - (endPosition - localPosition_);

            // Make sure to use a large enough vector to cross the negative
            // face. Bit overkill.
            const vector d(endPosition - localPosition_);
            const scalar magD(mag(d));
            if (magD > ROOTVSMALL)
            {
                // Get overall mesh bounding box
                treeBoundBox meshBb(mesh().bounds());
                // Extend to make 3D
                meshBb.inflate(ROOTSMALL);

                // Create vector guaranteed to cross mesh bounds
                projectedEndPosition = localPosition_ - meshBb.mag()*d/magD;

                // Clip to mesh bounds
                point intPt;
                direction intPtBits;
                bool ok = meshBb.intersects
                (
                    projectedEndPosition,
                    localPosition_ - projectedEndPosition,
                    projectedEndPosition,
                    localPosition_,
                    intPt,
                    intPtBits
                );
                if (ok)
                {
                    // Should always be the case
                    projectedEndPosition = intPt;
                }
            }
        }

        // Remove normal component
        {
            const point& basePt = mesh().points()[tri[0]];
            projectedEndPosition -= ((projectedEndPosition - basePt)&n)*n;
        }


        bool doTrack = false;
        if (meshEdgeStart_ == -1 && diagEdge_ == -1)
        {
            // We're starting and not yet on an edge.
            doTrack = true;
        }
        else
        {
            // See if the current triangle has got a point on the
            // correct side of the edge.
            doTrack = isTriAlongTrack(n, projectedEndPosition);
        }


        if (doTrack)
        {
            // Track across triangle. Return triangle edge crossed.
            label triEdgei = -1;
            trackFraction = trackFaceTri(n, projectedEndPosition, triEdgei);

            if (triEdgei == -1)
            {
                // Reached endpoint
                //checkInside();
                diagEdge_ = -1;
                meshEdgeStart_ = -1;
                return trackFraction;
            }

            const tetIndices ti(currentTetIndices());
            const triFace trif(ti.triIs(mesh(), false));
            // Triangle (faceTriIs) gets constructed from
            //    f[faceBasePtI_],
            //    f[facePtAI_],
            //    f[facePtBI_]
            //
            // So edge indices are:
            // 0 : edge between faceBasePtI_ and facePtAI_
            // 1 : edge between facePtAI_ and facePtBI_ (is always a real edge)
            // 2 : edge between facePtBI_ and faceBasePtI_

            const Foam::face& f = mesh().faces()[ti.face()];
            const label fp0 = trif[0];

            if (triEdgei == 0)
            {
                if (trif[1] == f.fcIndex(fp0))
                {
                    //Pout<< "Real edge." << endl;
                    diagEdge_ = -1;
                    meshEdgeStart_ = fp0;
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(cloud, td, trackFraction);
                }
                else if (trif[1] == f.rcIndex(fp0))
                {
                    //Note: should not happen since boundary face so owner
                    //Pout<< "Real edge." << endl;
                    FatalErrorInFunction
                        << abort(FatalError);

                    diagEdge_ = -1;
                    meshEdgeStart_ = f.rcIndex(fp0);
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(cloud, td, trackFraction);
                }
                else
                {
                    // Get index of triangle on other side of edge.
                    diagEdge_ = trif[1] - fp0;
                    if (diagEdge_ < 0)
                    {
                        diagEdge_ += f.size();
                    }
                    meshEdgeStart_ = -1;
                    //checkOnEdge();
                    crossDiagonalEdge();
                }
            }
            else if (triEdgei == 1)
            {
                //Pout<< "Real edge." << endl;
                diagEdge_ = -1;
                meshEdgeStart_ = trif[1];
                //checkOnEdge();
                crossEdgeConnectedFace(currentEdge());
                patchInteraction(cloud, td, trackFraction);
            }
            else // if (triEdgei == 2)
            {
                if (trif[2] == f.rcIndex(fp0))
                {
                    //Pout<< "Real edge." << endl;
                    diagEdge_ = -1;
                    meshEdgeStart_ = trif[2];
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(cloud, td, trackFraction);
                }
                else if (trif[2] == f.fcIndex(fp0))
                {
                    //Note: should not happen since boundary face so owner
                    //Pout<< "Real edge." << endl;
                    FatalErrorInFunction << abort(FatalError);

                    diagEdge_ = -1;
                    meshEdgeStart_ = fp0;
                    //checkOnEdge();
                    crossEdgeConnectedFace(currentEdge());
                    patchInteraction(cloud, td, trackFraction);
                }
                else
                {
                    //Pout<< "Triangle edge." << endl;
                    // Get index of triangle on other side of edge.
                    diagEdge_ = trif[2] - fp0;
                    if (diagEdge_ < 0)
                    {
                        diagEdge_ += f.size();
                    }
                    meshEdgeStart_ = -1;
                    //checkOnEdge();
                    crossDiagonalEdge();
                }
            }
        }
        else
        {
            // Current tet is not the right one. Check the neighbour tet.

            if (meshEdgeStart_ != -1)
            {
                // Particle is on mesh edge so change into other face on cell
                crossEdgeConnectedFace(currentEdge());
                //checkOnEdge();
                patchInteraction(cloud, td, trackFraction);
            }
            else
            {
                // Particle is on diagonal edge so change into the other
                // triangle.
                crossDiagonalEdge();
                //checkOnEdge();
            }
        }
    }

    //checkInside();

    return trackFraction;
}


template<class TrackCloudType>
void Foam::wallBoundedParticle::hitProcessorPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    // Switch particle
    td.switchProcessor = true;

    // Adapt edgeStart_ for other side.
    // E.g. if edgeStart_ is 1 then the edge is between vertex 1 and 2 so
    // on the other side between 2 and 3 so edgeStart_ should be
    // f.size()-edgeStart_-1.

    const Foam::face& f = mesh().faces()[face()];

    if (meshEdgeStart_ != -1)
    {
        meshEdgeStart_ = f.size() - meshEdgeStart_-1;
    }
    else
    {
        // diagEdge_ is relative to faceBasePt
        diagEdge_ = f.size() - diagEdge_;
    }
}


template<class TrackCloudType>
void Foam::wallBoundedParticle::hitWallPatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{}


template<class TrackCloudType>
void Foam::wallBoundedParticle::readFields(TrackCloudType& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<point> localPosition
    (
        c.newIOobject("position", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, localPosition);

    IOField<label> meshEdgeStart
    (
        c.newIOobject("meshEdgeStart", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, meshEdgeStart);

    IOField<label> diagEdge
    (
        c.newIOobject("diagEdge", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, diagEdge);

    label i = 0;
    for (wallBoundedParticle& p : c)
    {
        p.localPosition_ = localPosition[i];
        p.meshEdgeStart_ = meshEdgeStart[i];
        p.diagEdge_ = diagEdge[i];

        ++i;
    }
}


template<class TrackCloudType>
void Foam::wallBoundedParticle::writeFields(const TrackCloudType& c)
{
    particle::writeFields(c);

    const label np = c.size();

    IOField<point> localPosition
    (
        c.newIOobject("position", IOobject::NO_READ),
        np
    );
    IOField<label> meshEdgeStart
    (
        c.newIOobject("meshEdgeStart", IOobject::NO_READ),
        np
    );
    IOField<label> diagEdge
    (
        c.newIOobject("diagEdge", IOobject::NO_READ),
        np
    );

    label i = 0;
    for (const wallBoundedParticle& p : c)
    {
        localPosition[i] = p.localPosition_;
        meshEdgeStart[i] = p.meshEdgeStart_;
        diagEdge[i] = p.diagEdge_;

        ++i;
    }

    localPosition.write();
    meshEdgeStart.write();
    diagEdge.write();
}


// ************************************************************************* //

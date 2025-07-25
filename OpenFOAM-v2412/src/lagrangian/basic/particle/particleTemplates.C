/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2020 OpenFOAM Foundation
    Copyright (C) 2016-2020,2022 OpenCFD Ltd.
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

#include "IOPosition.H"

#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"
#include "cyclicACMIPolyPatch.H"
#include "processorPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wallPolyPatch.H"
#include "wedgePolyPatch.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::particle::writePropertyName
(
    Ostream& os,
    const word& name,
    const word& delim
)
{
    if (pTraits<Type>::nComponents == 1)
    {
        os  << name;
    }
    else
    {
        os  << '(';
        for (int i = 0; i < pTraits<Type>::nComponents; ++i)
        {
            if (i) os  << delim;

            os   << name << Foam::name(i);
        }
        os  << ')';
    }
}


template<class Type>
void Foam::particle::writeProperty
(
    Ostream& os,
    const word& name,
    const Type& value,
    const bool nameOnly,
    const word& delim,
    const wordRes& filters
)
{
    if (!filters.empty() && !filters.match(name))
    {
        return;
    }

    os  << delim;
    if (nameOnly)
    {
        writePropertyName<Type>(os, name, delim);
    }
    else
    {
        os  << value;
    }
}


template<class Type>
void Foam::particle::writeProperty
(
    Ostream& os,
    const word& name,
    const Field<Type>& values,
    const bool nameOnly,
    const word& delim,
    const wordRes& filters
)
{
    if (!filters.empty() && !filters.match(name))
    {
        return;
    }

    if (nameOnly)
    {
        os  << delim;
        os  << "N(";
        if (values.size())
        {
            forAll(values, i)
            {
                if (i) os  << delim;
                const word tag(name + Foam::name(i));
                writePropertyName<Type>(os, tag, delim);
            }
        }
        else
        {
            os  << name;
        }
        os << ')';
    }
    else
    {
        os  << delim << values;
    }
}


template<class TrackCloudType>
void Foam::particle::readFields(TrackCloudType& c)
{
    const bool readOnProc = c.size();

    IOobject procIO(c.newIOobject("origProcId", IOobject::MUST_READ));

    const bool haveFile = procIO.typeHeaderOk<IOField<label>>(true);

    IOField<label> origProcId(procIO, readOnProc && haveFile);
    c.checkFieldIOobject(c, origProcId);

    IOField<label> origId
    (
        c.newIOobject("origId", IOobject::MUST_READ),
        readOnProc && haveFile
    );
    c.checkFieldIOobject(c, origId);

    label i = 0;
    for (particle& p : c)
    {
        p.origProc_ = origProcId[i];
        p.origId_ = origId[i];

        ++i;
    }
}


template<class TrackCloudType>
void Foam::particle::writeFields(const TrackCloudType& c)
{
    const label np = c.size();
    const bool writeOnProc = c.size();

    if (writeLagrangianCoordinates)
    {
        IOPosition<TrackCloudType> ioP(c);
        ioP.write(writeOnProc);
    }
    else if (!writeLagrangianPositions)
    {
        FatalErrorInFunction
            << "Must select coordinates and/or positions" << nl
            << exit(FatalError);
    }

    // Optionally write positions file in v1706 format and earlier
    if (writeLagrangianPositions)
    {
        IOPosition<TrackCloudType> ioP
        (
            c,
            cloud::geometryType::POSITIONS
        );
        ioP.write(writeOnProc);
    }

    IOField<label> origProc
    (
        c.newIOobject("origProcId", IOobject::NO_READ),
        np
    );
    IOField<label> origId
    (
        c.newIOobject("origId", IOobject::NO_READ),
        np
    );

    label i = 0;
    for (const particle& p : c)
    {
        origProc[i] = p.origProc_;
        origId[i] = p.origId_;

        ++i;
    }

    origProc.write(writeOnProc);
    origId.write(writeOnProc);
}


template<class CloudType>
void Foam::particle::readObjects(CloudType& c, const objectRegistry& obr)
{
    typedef typename CloudType::parcelType parcelType;

    const auto* positionPtr = cloud::findIOPosition(obr);

    const label np = c.size();
    const label newNp = (positionPtr ? positionPtr->size() : 0);

    // Remove excess parcels
    for (label i = newNp; i < np; ++i)
    {
        parcelType* p = c.last();

        c.deleteParticle(*p);
    }

    if (newNp)
    {
        const auto& position = *positionPtr;

        const auto& origProcId = cloud::lookupIOField<label>("origProc", obr);
        const auto& origId = cloud::lookupIOField<label>("origId", obr);

        // Create new parcels
        for (label i = np; i < newNp; ++i)
        {
            c.addParticle(new parcelType(c.pMesh(), position[i], -1));
        }

        label i = 0;
        for (particle& p : c)
        {
            p.origProc_ = origProcId[i];
            p.origId_ = origId[i];

            if (i < np)
            {
                // Use relocate for old particles, not new ones
                p.relocate(position[i]);
            }

            ++i;
        }
    }
}


template<class CloudType>
void Foam::particle::writeObjects(const CloudType& c, objectRegistry& obr)
{
    const label np = c.size();

    auto& origProc = cloud::createIOField<label>("origProc", np, obr);
    auto& origId = cloud::createIOField<label>("origId", np, obr);
    auto& position = cloud::createIOField<point>("position", np, obr);

    label i = 0;
    for (const particle& p : c)
    {
        origProc[i] = p.origProc_;
        origId[i] = p.origId_;
        position[i] = p.position();

        ++i;
    }
}


template<class TrackCloudType>
void Foam::particle::hitBoundaryFace
(
    const vector& displacement,
    TrackCloudType& cloud,
    trackingData& td
)
{
    typename TrackCloudType::particleType& p =
        static_cast<typename TrackCloudType::particleType&>(*this);
    typename TrackCloudType::particleType::trackingData& ttd =
        static_cast<typename TrackCloudType::particleType::trackingData&>(td);

    if (!p.hitPatch(cloud, ttd))
    {
        const polyPatch& patch = mesh_.boundaryMesh()[p.patch()];

        if (isA<wedgePolyPatch>(patch))
        {
            p.hitWedgePatch(cloud, ttd);
        }
        else if (isA<symmetryPlanePolyPatch>(patch))
        {
            p.hitSymmetryPlanePatch(cloud, ttd);
        }
        else if (isA<symmetryPolyPatch>(patch))
        {
            p.hitSymmetryPatch(cloud, ttd);
        }
        else if (isA<cyclicPolyPatch>(patch))
        {
            p.hitCyclicPatch(cloud, ttd);
        }
        else if (isA<cyclicACMIPolyPatch>(patch))
        {
            p.hitCyclicACMIPatch(cloud, ttd, displacement);
        }
        else if (isA<cyclicAMIPolyPatch>(patch))
        {
            p.hitCyclicAMIPatch(cloud, ttd, displacement);
        }
        else if (isA<processorPolyPatch>(patch))
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


template<class TrackCloudType>
void Foam::particle::hitFace
(
    const vector& displacement,
    TrackCloudType& cloud,
    trackingData& td
)
{
    if (!onFace())
    {
        return;
    }
    else if (onInternalFace())
    {
        changeCell();
    }
    else if (onBoundaryFace())
    {
        // If on duplicate baffle (e.g. cyclicACMI) change to lowest numbered
        // patch ('masterpatch')
        changeToMasterPatch();

        hitBoundaryFace(displacement, cloud, td);
    }
}


template<class TrackCloudType>
void Foam::particle::trackToAndHitFace
(
    const vector& direction,
    const scalar fraction,
    TrackCloudType& cloud,
    trackingData& td
)
{
    trackToFace(direction, fraction);

    hitFace(direction, cloud, td);
}


template<class TrackCloudType>
bool Foam::particle::hitPatch(TrackCloudType&, trackingData&)
{
    return false;
}


template<class TrackCloudType>
void Foam::particle::hitWedgePatch(TrackCloudType& cloud, trackingData& td)
{
    FatalErrorInFunction
        << "Hitting a wedge patch should not be possible."
        << abort(FatalError);

    hitSymmetryPatch(cloud, td);
}


template<class TrackCloudType>
void Foam::particle::hitSymmetryPlanePatch
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    hitSymmetryPatch(cloud, td);
}


template<class TrackCloudType>
void Foam::particle::hitSymmetryPatch(TrackCloudType&, trackingData&)
{
    const vector nf = normal();

    transformProperties(I - 2.0*nf*nf);
}


template<class TrackCloudType>
void Foam::particle::hitCyclicPatch(TrackCloudType&, trackingData&)
{
    const cyclicPolyPatch& cpp =
        static_cast<const cyclicPolyPatch&>(mesh_.boundaryMesh()[patch()]);
    const cyclicPolyPatch& receiveCpp = cpp.neighbPatch();
    const label receiveFacei = receiveCpp.whichFace(facei_);

    // Set the topology
    facei_ = tetFacei_ = cpp.transformGlobalFace(facei_);
    celli_ = mesh_.faceOwner()[facei_];
    // See note in correctAfterParallelTransfer for tetPti addressing ...
    tetPti_ = mesh_.faces()[tetFacei_].size() - 1 - tetPti_;

    // Reflect to account for the change of triangle orientation in the new cell
    reflect();

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[receiveFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[receiveFacei]
        );
        transformProperties(-s);
    }
}


template<class TrackCloudType>
void Foam::particle::hitCyclicAMIPatch
(
    TrackCloudType&,
    trackingData& td,
    const vector& displacement
)
{
    vector pos = position();

    const cyclicAMIPolyPatch& cpp =
        static_cast<const cyclicAMIPolyPatch&>(mesh_.boundaryMesh()[patch()]);
    const cyclicAMIPolyPatch& receiveCpp = cpp.neighbPatch();
    const label sendFacei = cpp.whichFace(facei_);
    const label receiveFacei = cpp.pointFace(sendFacei, displacement, pos);

    if (receiveFacei < 0)
    {
        // If the patch face of the particle is not known assume that the
        // particle is lost and mark it to be deleted.
        td.keepParticle = false;
        WarningInFunction
            << "Particle lost across " << cyclicAMIPolyPatch::typeName
            << " patches " << cpp.name() << " and " << receiveCpp.name()
            << " at position " << pos << endl;
    }

    // Set the topology
    facei_ = tetFacei_ = receiveFacei + receiveCpp.start();

    // Locate the particle on the receiving side
    vector displacementT = displacement;
    cpp.reverseTransformDirection(displacementT, sendFacei);

    // NOTE: The ray used to find the hit location accross the AMI might not
    // be consistent in the displacement direction. Therefore a particle can
    // be looping accross AMI patches indefinitely. Advancing the particle
    // trajectory inside the cell is a possible solution.
    const vector dispDir = cpp.fraction()*displacementT;
    stepFraction_ += cpp.fraction();
    locate
    (
        pos + dispDir,
        &displacementT,
        mesh_.faceOwner()[facei_],
        false,
        "Particle crossed between " + cyclicAMIPolyPatch::typeName +
        " patches " + cpp.name() + " and " + receiveCpp.name() +
        " to a location outside of the mesh."
    );

    // The particle must remain associated with a face for the tracking to
    // register as incomplete
    facei_ = tetFacei_;

    // Transform the properties
    if (!receiveCpp.parallel())
    {
        const tensor& T =
        (
            receiveCpp.forwardT().size() == 1
          ? receiveCpp.forwardT()[0]
          : receiveCpp.forwardT()[receiveFacei]
        );
        transformProperties(T);
    }
    else if (receiveCpp.separated())
    {
        const vector& s =
        (
            (receiveCpp.separation().size() == 1)
          ? receiveCpp.separation()[0]
          : receiveCpp.separation()[receiveFacei]
        );
        transformProperties(-s);
    }

    //if (onBoundaryFace())
    {
//         vector receiveNormal, receiveDisplacement;
//         patchData(receiveNormal, receiveDisplacement);
//
//         if (((displacementT - fraction*receiveDisplacement)&receiveNormal) > 0)
//         {
//             td.keepParticle = false;
//             WarningInFunction
//                 << "Particle transfer from " << cyclicAMIPolyPatch::typeName
//                 << " patches " << cpp.name() << " to " << receiveCpp.name()
//                 << " failed at position " << pos << " and with displacement "
//                 << (displacementT - fraction*receiveDisplacement) << nl
//                 << "    The displacement points into both the source and "
//                 << "receiving faces, so the tracking cannot proceed" << nl
//                 << "    The particle has been removed" << nl << endl;
//             return;
//         }
    }
}


template<class TrackCloudType>
void Foam::particle::hitCyclicACMIPatch
(
    TrackCloudType& cloud,
    trackingData& td,
    const vector& displacement
)
{
    const cyclicACMIPolyPatch& cpp =
        static_cast<const cyclicACMIPolyPatch&>(mesh_.boundaryMesh()[patch()]);

    const label localFacei = cpp.whichFace(facei_);

    // If the mask is within the patch tolerance at either end, then we can
    // assume an interaction with the appropriate part of the ACMI pair.
    const scalar mask = cpp.mask()[localFacei];
    bool couple = mask >= 1 - cpp.tolerance();
    bool nonOverlap = mask <= cpp.tolerance();

    // If the mask is an intermediate value, then we search for a location on
    // the other side of the AMI. If we can't find a location, then we assume
    // that we have hit the non-overlap patch.
    if (!couple && !nonOverlap)
    {
        vector pos = position();
        couple = cpp.pointFace(localFacei, displacement, pos) >= 0;
        nonOverlap = !couple;
    }

    if (couple)
    {
        hitCyclicAMIPatch(cloud, td, displacement);
    }
    else
    {
        // Move to the face associated with the non-overlap patch and redo the
        // face interaction.
        tetFacei_ = facei_ = cpp.nonOverlapPatch().start() + localFacei;
        // Bypass all the changeToMasterPatch logic so as not to change
        // the patch back to cyclicACMI
        hitBoundaryFace(displacement, cloud, td);
    }
}


template<class TrackCloudType>
void Foam::particle::hitProcessorPatch(TrackCloudType&, trackingData&)
{}


template<class TrackCloudType>
void Foam::particle::hitWallPatch(TrackCloudType&, trackingData&)
{}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "polySurface.H"
#include "Time.H"
#include "ModifiableMeshedSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polySurface, 0);
}

const Foam::word Foam::polySurface::pointDataName("PointData");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polySurface::calculateZoneIds(const UList<surfZone>& zones)
{
    if (returnReduceAnd(zones.empty()))
    {
        zoneIds_.clear();
        return;
    }

    // Extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently

    zoneIds_.resize_nocopy(size());

    label off = 0;
    for (const surfZone& zn : zones)
    {
        SubList<label>(zoneIds_, zn.size(), off) = zn.index();
        off += zn.size();
    }

    if (off < size())
    {
        WarningInFunction
            << "More faces " << size() << " than zones " << off << endl;

        zoneIds_.slice(off) = (zones.empty() ? 0 : zones.back().index());
    }
    else if (size() < off)
    {
        FatalErrorInFunction
            << "More zones " << off << " than faces " << size()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polySurface::polySurface(const IOobject& io, bool doCheckIn)
:
    objectRegistry
    (
        IOobject
        (
            io.name(),
            io.db().time().constant(),
            io.db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::REGISTER
        )
    ),
    MeshReference(faceList(), pointField())
{
    // Created without a point field sub-registry

    if (doCheckIn)
    {
        this->store();
    }
}


Foam::polySurface::polySurface
(
    const word& surfName,
    const objectRegistry& obr,
    bool doCheckIn
)
:
    polySurface
    (
        IOobject
        (
            surfName,
            obr.time().constant(),
            obr,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::REGISTER
        ),
        doCheckIn
    )
{}


Foam::polySurface::polySurface
(
    const IOobject& io,
    const MeshedSurface<face>& surf,
    bool doCheckIn
)
:
    polySurface(io, doCheckIn)
{
    copySurface(surf);
}


Foam::polySurface::polySurface
(
    const IOobject& io,
    MeshedSurface<face>&& surf,
    bool doCheckIn
)
:
    polySurface(io, doCheckIn)
{
    transfer(surf);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::polySurface& Foam::polySurface::New
(
    const word& surfName,
    const objectRegistry& obr
)
{
    auto* ptr = obr.getObjectPtr<polySurface>(surfName);

    if (!ptr)
    {
        ptr = new polySurface(surfName, obr);
        regIOobject::store(ptr);
    }

    return *ptr;
}


bool Foam::polySurface::Delete
(
    const word& surfName,
    const objectRegistry& obr
)
{
    return obr.checkOut(obr.getObjectPtr<polySurface>(surfName));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polySurface::~polySurface()
{
    ///  clearOut(); // Clear addressing
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::polySurface::nFaceData() const
{
    // Do not include sub-registry in the count
    if (objectRegistry::foundObject<objectRegistry>(pointDataName))
    {
        return (objectRegistry::size() - 1);
    }

    return objectRegistry::size();
}


Foam::label Foam::polySurface::nPointData() const
{
    const auto* subreg =
        objectRegistry::cfindObject<objectRegistry>(pointDataName);

    if (subreg)
    {
        return subreg->size();
    }

    return 0;
}


const Foam::objectRegistry& Foam::polySurface::faceData() const
{
    return static_cast<const objectRegistry&>(*this);
}


const Foam::objectRegistry& Foam::polySurface::pointData() const
{
    // Force create on access
    return objectRegistry::subRegistry(pointDataName, true);
}


Foam::polySurface::FieldAssociation
Foam::polySurface::queryFieldAssociation(const word& fieldName) const
{
    unsigned where(FieldAssociation::NO_DATA);

    const objectRegistry* subreg = nullptr;

    // Face Data
    {
        subreg = this;

        if (subreg && subreg->contains(fieldName))
        {
            where |= FieldAssociation::FACE_DATA;
        }
    }

    // Point Data
    {
        subreg = cfindObject<objectRegistry>(pointDataName);

        if (subreg && subreg->contains(fieldName))
        {
            where |= FieldAssociation::POINT_DATA;
        }
    }

    return FieldAssociation(where);
}


const Foam::regIOobject* Foam::polySurface::findFieldObject
(
    const word& fieldName,
    enum FieldAssociation association
) const
{
    const unsigned where(association);


    const regIOobject* ioptr = nullptr;
    const objectRegistry* subreg = nullptr;

    // Face Data
    if (where & FieldAssociation::FACE_DATA)
    {
        subreg = this;

        if (subreg && (ioptr = subreg->cfindObject<regIOobject>(fieldName)))
        {
            return ioptr;
        }
    }

    // Point Data
    if (where & FieldAssociation::POINT_DATA)
    {
        subreg = cfindObject<objectRegistry>(pointDataName);

        if (subreg && (ioptr = subreg->cfindObject<regIOobject>(fieldName)))
        {
            return ioptr;
        }
    }

    return ioptr;
}


void Foam::polySurface::copySurface
(
    const pointField& points,
    const faceList& faces,
    bool unused
)
{
    clearOut(); // Clear addressing

    if
    (
        this->nPoints() != points.size()
     || this->nFaces() != faces.size()
    )
    {
        // Geometry changed
        clearFields();
    }

    this->storedPoints() = points;
    this->storedFaces() = faces;

    zoneIds_.clear();

    // if (validate)
    // {
    //     checkZones();
    // }
}


void Foam::polySurface::copySurface
(
    const meshedSurf& surf,
    bool unused
)
{
    clearOut(); // Clear addressing

    if
    (
        this->nPoints() != surf.points().size()
     || this->nFaces() != surf.faces().size()
    )
    {
        // Geometry changed
        clearFields();
    }

    this->storedPoints() = surf.points();
    this->storedFaces() = surf.faces();

    zoneIds_ = surf.zoneIds();

    // if (validate)
    // {
    //     checkZones();
    // }
}


void Foam::polySurface::copySurface
(
    const MeshedSurface<face>& surf,
    bool unused
)
{
    clearOut(); // Clear addressing

    if
    (
        this->nPoints() != surf.points().size()
     || this->nFaces() != surf.surfFaces().size()
    )
    {
        // Geometry changed
        clearFields();
    }

    this->storedPoints() = surf.points();
    this->storedFaces() = surf.surfFaces();

    calculateZoneIds(surf.surfZones());

    // if (validate)
    // {
    //     checkZones();
    // }
}


void Foam::polySurface::transfer
(
    pointField&& points,
    faceList&& faces,
    labelList&& zoneIds
)
{
    clearOut(); // Clear addressing
    clearFields();

    this->storedPoints().transfer(points);
    this->storedFaces().transfer(faces);
    zoneIds_.transfer(zoneIds);
}


void Foam::polySurface::transfer
(
    MeshedSurface<face>& surf,
    bool unused
)
{
    clearOut(); // Clear addressing
    clearFields();

    ModifiableMeshedSurface<face> tsurf(std::move(surf));

    this->storedPoints().transfer(tsurf.storedPoints());
    this->storedFaces().transfer(tsurf.storedFaces());

    calculateZoneIds(tsurf.surfZones());

    // if (validate)
    // {
    //     checkZones();
    // }
}


// void Foam::polySurface::checkZones()
// {
//     // Extra safety, ensure we have at some zones
//     // and they cover all the faces - fix start silently
//
//     if (surfZones_.size() <= 1)
//     {
//         removeZones();
//         return;
//     }
//
//     label count = 0;
//     for (surfZone& zn : surfZones_)
//     {
//         zn.start() = count;
//         count += zn.size();
//     }
//
//     if (count < nFaces())
//     {
//         WarningInFunction
//             << "More faces " << nFaces() << " than zones " << count
//             << " ... extending final zone"
//             << endl;
//
//         surfZones_.back().size() += count - nFaces();
//     }
//     else if (size() < count)
//     {
//         FatalErrorInFunction
//             << "More zones " << count << " than faces " << nFaces()
//             << exit(FatalError);
//     }
// }


// * * * * * * * * * * * * * * * Specializations * * * * * * * * * * * * * * //

namespace Foam
{

template<>
const regIOobject* polySurface::findFieldObject<polySurfaceGeoMesh>
(
    const word& fieldName
) const
{
    // Face Data (main registry)
    return cfindObject<regIOobject>(fieldName);
}


template<>
const regIOobject* polySurface::findFieldObject<polySurfacePointGeoMesh>
(
    const word& fieldName
) const
{
    // Point Data (sub-registry)

    const auto* subreg =
        objectRegistry::cfindObject<objectRegistry>(pointDataName);

    if (subreg)
    {
        return subreg->cfindObject<regIOobject>(fieldName);
    }

    return nullptr;
}



template<>
const objectRegistry* polySurface::whichRegistry<polySurfaceGeoMesh>
(
    const word& fieldName
) const
{
    // Face Data (main registry)
    const objectRegistry* subreg = this;

    if (subreg->contains(fieldName))
    {
        return subreg;
    }

    return nullptr;
}


template<>
const objectRegistry* polySurface::whichRegistry<polySurfacePointGeoMesh>
(
    const word& fieldName
) const
{
    // Point Data (sub registry)

    const auto* subreg =
        objectRegistry::cfindObject<objectRegistry>(pointDataName);

    if (subreg && subreg->contains(fieldName))
    {
        return subreg;
    }

    return nullptr;
}

} // End of namespace

// ************************************************************************* //

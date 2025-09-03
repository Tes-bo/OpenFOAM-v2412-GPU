/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "surfMesh.H"
#include "meshedSurf.H"
#include "MeshedSurfaceProxy.H"

#include "Time.H"
#include "OSspecific.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMesh, 0);
}

Foam::word Foam::surfMesh::meshSubDir = "surfMesh";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// void Foam::surfMesh::oneZone()
// {
//     word zoneName;
//
//     if (surfZones_.size())
//     {
//         zoneName = surfZones_[0].name();
//     }
//     if (zoneName.empty())
//     {
//         zoneName = "zone0";
//     }
//
//     // Set single default zone with nFaces
//     surfZones_.resize(1);
//     surfZones_[0] = surfZone(zoneName, nFaces());
// }


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::surfMesh::updateRefs()
{
    // Synchronize UList reference to the faces
    static_cast<MeshReference&>(*this).shallowCopy
    (
        this->storedFaces()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMesh::surfMesh(const IOobject& io)
:
    surfMesh(io, io.name())
{}


Foam::surfMesh::surfMesh(const IOobject& io, const word& surfName)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE
        ),
        IOobject
        (
            "faces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobjectOption::MUST_READ,
            IOobjectOption::NO_WRITE
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            time().findInstance(meshDir(), "surfZones"),
            meshSubDir,
            *this,
            IOobjectOption::LAZY_READ,
            IOobjectOption::NO_WRITE
        )
    )
{}


Foam::surfMesh::surfMesh
(
    const word& surfName,
    const objectRegistry& obr
)
:
    surfaceRegistry(obr, surfName),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE
        ),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE
        )
    )
{}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    const MeshedSurface<face>& surf,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobjectOption::NO_READ,
            io.writeOpt()
        ),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobjectOption::NO_READ,
            io.writeOpt()
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            IOobjectOption::NO_READ,
            io.writeOpt()
        )
    )
{
    DebugInfo
        <<"IOobject: " << io.path() << nl
        <<"  name: " << io.name()
        <<"  instance: " << io.instance()
        <<"  local: " << io.local()
        <<"  dbDir: " << io.db().dbDir() << nl
        <<"creating surfMesh at instance " << instance() << endl;

    copySurface(surf);
}


Foam::surfMesh::surfMesh
(
    const IOobject& io,
    MeshedSurface<face>&& surf,
    const word& surfName
)
:
    surfaceRegistry(io.db(), (surfName.size() ? surfName : io.name())),
    Allocator
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        ),
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    ),
    MeshReference(this->storedIOFaces(), this->storedIOPoints()),

    surfZones_
    (
        IOobject
        (
            "surfZones",
            instance(),
            meshSubDir,
            *this,
            io.readOpt(),
            io.writeOpt()
        )
    )
{
    DebugInfo
        <<"IOobject: " << io.path() << nl
        <<" name: " << io.name()
        <<" instance: " << io.instance()
        <<" local: " << io.local()
        <<" dbDir: " << io.db().dbDir() << nl
        <<"creating surfMesh at instance " << instance() << nl
        <<"timeName: " << instance() << endl;

    transfer(surf);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMesh::~surfMesh()
{
    clearOut(); // Clear addressing
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::copySurface
(
    const pointField& points,
    const faceList& faces,
    bool validate
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

    this->storedIOPoints() = points;
    this->storedIOFaces() = faces;
    surfZones_.clear();

    this->updateRefs();

    // No zones
}


void Foam::surfMesh::copySurface
(
    const meshedSurf& surf,
    bool validate
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

    this->storedIOPoints() = surf.points();
    this->storedIOFaces() = surf.faces();
    surfZones_.clear();

    this->updateRefs();

    // No zones
}


void Foam::surfMesh::copySurface
(
    const MeshedSurface<face>& surf,
    bool validate
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

    this->storedIOPoints() = surf.points();
    this->storedIOFaces() = surf.surfFaces();
    surfZones_ = surf.surfZones();

    this->updateRefs();

    if (validate)
    {
        checkZones();
    }
}


void Foam::surfMesh::transfer
(
    MeshedSurface<face>& surf,
    bool validate
)
{
    clearOut(); // Clear addressing

    this->storedPoints().transfer(surf.storedPoints());
    this->storedFaces().transfer(surf.storedFaces());
    this->storedZones().transfer(surf.storedZones());

    this->updateRefs();

    if (validate)
    {
        checkZones();
    }
}


Foam::autoPtr<Foam::MeshedSurface<Foam::face>>
Foam::surfMesh::releaseGeom()
{
    clearOut(); // Clear addressing
    clearFields();

    // Start with an empty geometry
    auto aptr = autoPtr<MeshedSurface<face>>::New();

    // Transfer in content
    aptr->storedPoints().transfer(this->storedPoints());
    aptr->storedFaces().transfer(this->storedFaces());
    aptr->storedZones().transfer(this->storedZones());

    this->updateRefs(); // This may not be needed...

    return aptr;
}


Foam::fileName Foam::surfMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const Foam::fileName& Foam::surfMesh::pointsInstance() const
{
    return this->storedIOPoints().instance();
}


const Foam::fileName& Foam::surfMesh::facesInstance() const
{
    return this->storedIOFaces().instance();
}


Foam::label Foam::surfMesh::nPoints() const
{
    return this->points().size();
}


Foam::label Foam::surfMesh::nFaces() const
{
    return this->faces().size();
}


const Foam::pointField& Foam::surfMesh::points() const
{
    return this->storedIOPoints();
}


const Foam::faceList& Foam::surfMesh::faces() const
{
    return this->storedIOFaces();
}


void Foam::surfMesh::checkZones(const bool verbose)
{
    auto& zones = this->storedZones();

    if (zones.size() <= 1)
    {
        removeZones();
        return;
    }

    // Check that zones cover all faces exactly
    // - fix each start silently

    bool zonesTooBig(false);

    const label maxCount = this->nFaces();

    label start = 0;
    for (surfZone& zn : zones)
    {
        zn.start() = start;
        start += zn.size();
        if (start > maxCount)
        {
            zonesTooBig = true;  // Zones exceed what is required
            zn.size() = (maxCount - zn.start());
            start = (zn.start() + zn.size());
        }
    }

    if (!zones.empty())
    {
        surfZone& zn = zones.back();

        if ((zn.start() + zn.size()) < maxCount)
        {
            // Zones address less than expected - extend final zone
            zn.size() += maxCount - zn.start();

            if (verbose)
            {
                WarningInFunction
                    << "Surface has more faces " << maxCount
                    << " than zone addressing ... extending final zone" << nl;
            }
        }
        else if (zonesTooBig)
        {
            if (verbose)
            {
                WarningInFunction
                    << "Surface has more zone addressing than faces "
                    << maxCount
                    << " ... trucated/resized accordingly" << nl;
            }
        }
    }
}


// Add boundary patches. Constructor helper
void Foam::surfMesh::addZones
(
    const surfZoneList& zones,
    const bool validate
)
{
    removeZones();

    forAll(surfZones_, zonei)
    {
        surfZones_[zonei] = surfZone(zones[zonei], zonei);
    }

    if (validate)
    {
        checkZones();
    }
}


// Remove all files and some subdirs (eg, sets)
void Foam::surfMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = db().path()/instanceDir/meshSubDir;

    rm(meshFilesPath/"points");
    rm(meshFilesPath/"faces");
    rm(meshFilesPath/"surfZones");
}


void Foam::surfMesh::removeFiles() const
{
    removeFiles(instance());
}


void Foam::surfMesh::write
(
    const fileName& name,
    IOstreamOption streamOpt,
    const dictionary& options
) const
{
    write(name, name.ext(), streamOpt, options);
}


void Foam::surfMesh::write
(
    const fileName& name,
    const word& fileType,
    IOstreamOption streamOpt,
    const dictionary& options
) const
{
    MeshedSurfaceProxy<face>
    (
        this->points(),
        this->faces(),
        this->surfZones()
    ).write(name, fileType, streamOpt, options);
}


// ************************************************************************* //

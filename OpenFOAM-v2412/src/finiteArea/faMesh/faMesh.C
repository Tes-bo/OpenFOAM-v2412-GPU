/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "faMesh.H"
#include "faMeshBoundaryHalo.H"
#include "faGlobalMeshData.H"
#include "Time.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "IndirectList.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "faMeshLduAddressing.H"
#include "processorFaPatch.H"
#include "wedgeFaPatch.H"
#include "faPatchData.H"
#include "registerSwitch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMesh, 0);

    int faMesh::geometryOrder_
    (
        debug::optimisationSwitch("fa:geometryOrder", 2)
    );
    registerOptSwitch
    (
        "fa:geometryOrder",
        int,
        faMesh::geometryOrder_
    );
}


const Foam::word Foam::faMesh::prefix_("finite-area");

Foam::word Foam::faMesh::meshSubDir("faMesh");

const int Foam::faMesh::quadricsFit_ = 0;  // Tuning (experimental)


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::word& Foam::faMesh::prefix() noexcept
{
    return prefix_;
}


Foam::fileName Foam::faMesh::dbDir(const word& areaRegion)
{
    if (areaRegion.empty() || areaRegion == polyMesh::defaultRegion)
    {
        return faMesh::prefix();
    }

    return (faMesh::prefix() / areaRegion);
}


Foam::fileName Foam::faMesh::dbDir
(
    const word& volRegion,
    const word& areaRegion
)
{
    return
    (
        polyMesh::regionName(volRegion)
      / faMesh::prefix()
      / polyMesh::regionName(areaRegion)
    );
}


Foam::fileName Foam::faMesh::dbDir
(
    const polyMesh& pMesh,
    const word& areaRegion
)
{
    return faMesh::dbDir(pMesh.regionName(), areaRegion);
}


Foam::fileName Foam::faMesh::meshDir(const word& areaRegion)
{
    if (areaRegion.empty() || areaRegion == polyMesh::defaultRegion)
    {
        return faMesh::meshSubDir;
    }

    return (areaRegion / faMesh::meshSubDir);
}


Foam::fileName Foam::faMesh::meshDir
(
    const word& volRegion,
    const word& areaRegion
)
{
    return
    (
        polyMesh::regionName(volRegion)
      / faMesh::prefix()
      / polyMesh::regionName(areaRegion)
      / faMesh::meshSubDir
    );
}


Foam::fileName Foam::faMesh::meshDir
(
    const polyMesh& pMesh,
    const word& areaRegion
)
{
    return faMesh::meshDir(pMesh.regionName(), areaRegion);
}


const Foam::objectRegistry* Foam::faMesh::registry(const polyMesh& pMesh)
{
    return pMesh.cfindObject<objectRegistry>(faMesh::prefix());
}


// const Foam::objectRegistry* Foam::faMesh::registry(const objectRegistry& obr)
// {
//     return obr.cfindObject<objectRegistry>(faMesh::prefix());
// }


const Foam::faMesh& Foam::faMesh::mesh
(
    const polyMesh& pMesh
)
{
    return faMesh::mesh(pMesh, polyMesh::defaultRegion);
}


const Foam::faMesh& Foam::faMesh::mesh
(
    const polyMesh& pMesh,
    const word& areaRegion
)
{
    const objectRegistry* obr = faMesh::registry(pMesh);

    if (!obr)
    {
        // Fallback - not really valid, but will fail at the next stage
        obr = &(pMesh.thisDb());
    }

    if (areaRegion.empty())
    {
        return obr->lookupObject<faMesh>(polyMesh::defaultRegion);
    }

    return obr->lookupObject<faMesh>(areaRegion);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Convert patch names to face labels. Preserve patch order
static labelList selectPatchFaces
(
    const polyBoundaryMesh& pbm,
    const wordRes& polyPatchNames
)
{
    const labelList patchIDs
    (
        pbm.indices(polyPatchNames, true)  // useGroups
    );

    if (patchIDs.empty())
    {
        FatalErrorInFunction
            << "No matching patches: " << polyPatchNames << nl
            << exit(FatalError);
    }

    label nFaceLabels = 0;
    for (const label patchi : patchIDs)
    {
        nFaceLabels += pbm[patchi].size();
    }

    labelList faceLabels(nFaceLabels);

    nFaceLabels = 0;
    for (const label patchi : patchIDs)
    {
        for (const label facei : pbm[patchi].range())
        {
            faceLabels[nFaceLabels] = facei;
            ++nFaceLabels;
        }
    }

    return faceLabels;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMesh::checkBoundaryEdgeLabelRange
(
    const labelUList& edgeLabels
) const
{
    label nErrors = 0;

    for (const label edgei : edgeLabels)
    {
        if (edgei < nInternalEdges_ || edgei >= nEdges_)
        {
            if (!nErrors++)
            {
                FatalErrorInFunction
                    << "Boundary edge label out of range "
                    << nInternalEdges_ << ".." << (nEdges_-1) << nl
                    << "   ";
            }

            FatalError<< ' ' << edgei;
        }
    }

    if (nErrors)
    {
        FatalError << nl << exit(FatalError);
    }
}


void Foam::faMesh::initPatch() const
{
    patchPtr_ = std::make_unique<uindirectPrimitivePatch>
    (
        UIndirectList<face>(mesh().faces(), faceLabels_),
        mesh().points()
    );

    // Could set some basic primitive data here...
    // nEdges_ = patchPtr_->nEdges();
    // nInternalEdges_ = patchPtr_->nInternalEdges();
    // nFaces_ = patchPtr_->size();
    // nPoints_ = patchPtr_->nPoints();
    bndConnectPtr_.reset(nullptr);
    haloMapPtr_.reset(nullptr);
    haloFaceCentresPtr_.reset(nullptr);
    haloFaceNormalsPtr_.reset(nullptr);
}


void Foam::faMesh::setPrimitiveMeshData()
{
    DebugInFunction << "Setting primitive data" << endl;

    const uindirectPrimitivePatch& bp = patch();
    const labelListList& edgeFaces = bp.edgeFaces();

    // Dimensions

    nEdges_ = bp.nEdges();
    nInternalEdges_ = bp.nInternalEdges();
    nFaces_ = bp.size();
    nPoints_ = bp.nPoints();

    edges_.resize(nEdges_);
    edgeOwner_.resize(nEdges_);
    edgeNeighbour_.resize(nInternalEdges_);

    // Internal edges
    for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
    {
        edges_[edgei] = bp.edges()[edgei];

        edgeOwner_[edgei] = edgeFaces[edgei][0];

        edgeNeighbour_[edgei] = edgeFaces[edgei][1];
    }

    // Continue with boundary edges
    label edgei = nInternalEdges_;

    for (const faPatch& p : boundary())
    {
        for (const label patchEdgei : p.edgeLabels())
        {
            edges_[edgei] = bp.edges()[patchEdgei];

            edgeOwner_[edgei] = edgeFaces[patchEdgei][0];

            ++edgei;
        }
    }
}


void Foam::faMesh::clearHalo() const
{
    DebugInFunction << "Clearing halo information" << endl;

    haloMapPtr_.reset(nullptr);
    haloFaceCentresPtr_.reset(nullptr);
    haloFaceNormalsPtr_.reset(nullptr);
}


void Foam::faMesh::clearGeomNotAreas() const
{
    DebugInFunction << "Clearing geometry" << endl;

    clearHalo();
    patchPtr_.reset(nullptr);
    polyPatchFacesPtr_.reset(nullptr);
    polyPatchIdsPtr_.reset(nullptr);
    bndConnectPtr_.reset(nullptr);
    SPtr_.reset(nullptr);
    patchStartsPtr_.reset(nullptr);
    LePtr_.reset(nullptr);
    magLePtr_.reset(nullptr);
    faceCentresPtr_.reset(nullptr);
    edgeCentresPtr_.reset(nullptr);
    faceAreaNormalsPtr_.reset(nullptr);
    edgeAreaNormalsPtr_.reset(nullptr);
    pointAreaNormalsPtr_.reset(nullptr);
    faceCurvaturesPtr_.reset(nullptr);
    edgeTransformTensorsPtr_.reset(nullptr);
}


void Foam::faMesh::clearGeom() const
{
    DebugInFunction << "Clearing geometry" << endl;

    clearGeomNotAreas();
    S0Ptr_.reset(nullptr);
    S00Ptr_.reset(nullptr);
    correctPatchPointNormalsPtr_.reset(nullptr);
}


void Foam::faMesh::clearAddressing() const
{
    DebugInFunction << "Clearing addressing" << endl;

    lduPtr_.reset(nullptr);
}


void Foam::faMesh::clearOut() const
{
    clearGeom();
    clearAddressing();
    globalMeshDataPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::faMesh::syncGeom()
{
    if (UPstream::parRun())
    {
        // areaCentres()
        if (faceCentresPtr_)
        {
            faceCentresPtr_->boundaryFieldRef()
                .evaluateCoupled<processorFaPatch>();
        }

        // faceAreaNormals()
        if (faceAreaNormalsPtr_)
        {
            faceAreaNormalsPtr_->boundaryFieldRef()
                .evaluateCoupled<processorFaPatch>();
        }
    }
}


bool Foam::faMesh::init(const bool doInit)
{
    if (doInit)
    {
        setPrimitiveMeshData();
    }

    // Create global mesh data
    if (UPstream::parRun())
    {
        (void)globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    syncGeom();

    return false;
}


// * * * * * * * * * * * * * Forwarding Constructors * * * * * * * * * * * * //

Foam::faMesh::faMesh
(
    const word& meshName,
    const polyMesh& pMesh,
    Foam::zero
)
:
    faMesh(meshName, pMesh, labelList())
{}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    Foam::zero
)
:
    faMesh(polyMesh::defaultRegion, pMesh, labelList())
{}


Foam::faMesh::faMesh(const polyMesh& pMesh, const bool doInit)
:
    faMesh(polyMesh::defaultRegion, pMesh, doInit)
{}


Foam::faMesh::faMesh
(
    const word& meshName,
    const faMesh& baseMesh,
    Foam::zero
)
:
    faMesh(meshName, baseMesh, labelList())
{}


Foam::faMesh::faMesh
(
    const faMesh& baseMesh,
    Foam::zero
)
:
    faMesh(polyMesh::defaultRegion, baseMesh, labelList())
{}


Foam::faMesh::faMesh
(
    const word& meshName,
    const faMesh& baseMesh,
    labelList&& faceLabels
)
:
    faMesh
    (
        meshName,
        baseMesh,
        std::move(faceLabels),
        static_cast<IOobjectOption>(baseMesh.thisDb())
    )
{}


Foam::faMesh::faMesh
(
    const faMesh& baseMesh,
    labelList&& faceLabels
)
:
    faMesh
    (
        polyMesh::defaultRegion,
        baseMesh,
        std::move(faceLabels),
        static_cast<IOobjectOption>(baseMesh.thisDb())
    )
{}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    labelList&& faceLabels,
    IOobjectOption ioOpt
)
:
    faMesh
    (
        polyMesh::defaultRegion,
        pMesh,
        std::move(faceLabels),
        ioOpt
    )
{}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    labelList&& faceLabels
)
:
    faMesh(polyMesh::defaultRegion, pMesh, std::move(faceLabels))
{}


Foam::faMesh::faMesh
(
    const polyPatch& pp,
    const bool doInit
)
:
    faMesh(polyMesh::defaultRegion, pp, doInit)
{}


Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const dictionary& faMeshDefinition,
    const bool doInit
)
:
    faMesh
    (
        polyMesh::defaultRegion,
        pMesh,
        faMeshDefinition,
        doInit
    )
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMesh::faMesh
(
    const word& meshName,
    const polyMesh& pMesh,
    const bool doInit
)
:
    faMeshRegistry(meshName, pMesh),
    faSchemes
    (
        faMesh::thisDb(),
        IOobjectOption::MUST_READ
    ),
    faSolution
    (
        faMesh::thisDb(),
        IOobjectOption::MUST_READ
    ),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            time().findInstance(meshDir(), "faceLabels"),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            // Allow boundary file that is newer than faceLabels
            time().findInstance
            (
                meshDir(),
                "faBoundary",
                IOobject::MUST_READ,
                faceLabels_.instance()
            ),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    comm_(UPstream::worldComm),
    curTimeIndex_(time().timeIndex())
{
    DebugInFunction << "Creating from IOobject" << endl;

    setPrimitiveMeshData();

    if (doInit)
    {
        faMesh::init(false);  // do not init lower levels
    }

    if (doInit)
    {
        // Read some optional fields
        // - logic as per fvMesh

        IOobject rio
        (
            "any-name",
            time().timeName(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::LAZY_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // Read old surface areas (if present)
        rio.resetHeader("S0");
        if (returnReduceOr(rio.typeHeaderOk<regIOobject>(false)))
        {
            S0Ptr_ = std::make_unique<DimensionedField<scalar, areaMesh>>
            (
                rio,
                *this,
                dimensionedScalar(dimArea, Zero)
            );
        }
    }
}


Foam::faMesh::faMesh
(
    const word& meshName,
    const polyMesh& pMesh,
    labelList&& faceLabels
)
:
    faMeshRegistry(meshName, pMesh),
    faSchemes
    (
        faMesh::thisDb(),
        IOobjectOption::MUST_READ
    ),
    faSolution
    (
        faMesh::thisDb(),
        IOobjectOption::MUST_READ
    ),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            pMesh.facesInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        std::move(faceLabels)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            faceLabels_.instance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    comm_(UPstream::worldComm),
    curTimeIndex_(time().timeIndex())
{
    // Not yet much for primitive mesh data possible...
    nPoints_ = 0;
    nEdges_ = 0;
    nInternalEdges_ = 0;
    nFaces_ = faceLabels_.size();

    // TDB: can we make a NO_READ readOption persistent for
    // faSchemes/faSolution? Or not needed anymore?
}

Foam::faMesh::faMesh
(
    const word& meshName,
    const polyMesh& pMesh,
    labelList&& faceLabels,
    IOobjectOption ioOpt
)
:
    faMeshRegistry(meshName, pMesh),
    faSchemes
    (
        faMesh::thisDb(),
        ioOpt.readOpt()
    ),
    faSolution
    (
        faMesh::thisDb(),
        ioOpt.readOpt()
    ),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            pMesh.facesInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        std::move(faceLabels)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            faceLabels_.instance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    comm_(UPstream::worldComm),
    curTimeIndex_(time().timeIndex())
{
    // Not yet much for primitive mesh data possible...
    nPoints_ = 0;
    nEdges_ = 0;
    nInternalEdges_ = 0;
    nFaces_ = faceLabels_.size();

    // TDB: can we make a NO_READ readOption persistent for
    // faSchemes/faSolution? Or not needed anymore?
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMesh::faMesh
(
    const word& meshName,
    const faMesh& baseMesh,
    labelList&& faceLabels,
    IOobjectOption ioOpt
)
:
    faMeshRegistry(meshName, baseMesh.mesh()),
    faSchemes
    (
        faMesh::thisDb(),
        ioOpt.readOpt(),
        static_cast<const dictionary*>(baseMesh.hasSchemes())
    ),
    faSolution
    (
        faMesh::thisDb(),
        ioOpt.readOpt(),
        static_cast<const dictionary*>(baseMesh.hasSolution())
    ),
    edgeInterpolation(*this),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            // Topological instance from polyMesh
            baseMesh.mesh().facesInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        std::move(faceLabels)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            faceLabels_.instance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        Foam::zero{}
    ),
    comm_(UPstream::worldComm),
    curTimeIndex_(time().timeIndex())
{
    // Not yet much for primitive mesh data possible...
    nPoints_ = 0;
    nEdges_ = 0;
    nInternalEdges_ = 0;
    nFaces_ = faceLabels_.size();
}


Foam::faMesh::faMesh
(
    const word& meshName,
    const polyPatch& pp,
    const bool doInit
)
:
    faMesh
    (
        meshName,
        pp.boundaryMesh().mesh(),
        identity(pp.range())
    )
{
    DebugInFunction << "Creating from polyPatch:" << pp.name() << endl;

    // Add single faPatch "default", but with processor connections
    faPatchList newPatches
    (
        createOnePatch("default")
    );

    addFaPatches(newPatches);

    setPrimitiveMeshData();

    if (doInit)
    {
        faMesh::init(false);  // do not init lower levels
    }
}


Foam::faMesh::faMesh
(
    const word& meshName,
    const polyMesh& pMesh,
    const dictionary& faMeshDefinition,
    const bool doInit
)
:
    faMesh
    (
        meshName,
        pMesh,
        selectPatchFaces
        (
            pMesh.boundaryMesh(),
            faMeshDefinition.get<wordRes>("polyMeshPatches")
        )
    )
{
    DebugInFunction << "Creating from definition (dictionary)" << endl;

    faPatchList newPatches
    (
        createPatchList
        (
            faMeshDefinition.subDict("boundary"),

            // Optional 'empty' patch
            faMeshDefinition.getOrDefault<word>("emptyPatch", word::null),

            // Optional specification for default patch
            faMeshDefinition.findDict("defaultPatch")
        )
    );

    addFaPatches(newPatches);

    if (doInit)
    {
        faMesh::init(false);  // do not init lower levels
    }

    if (doInit)
    {
        // Read old surface areas (if present)
        // - logic as per fvMesh

        IOobject rio
        (
            "any-name",
            time().timeName(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::LAZY_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        // Read old surface areas (if present)
        rio.resetHeader("S0");
        if (returnReduceOr(rio.typeHeaderOk<regIOobject>(false)))
        {
            S0Ptr_ = std::make_unique<DimensionedField<scalar, areaMesh>>
            (
                rio,
                *this,
                dimensionedScalar(dimArea, Zero)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faMesh::~faMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faSchemes* Foam::faMesh::hasSchemes() const
{
    return static_cast<const faSchemes*>(this);
}


const Foam::faSolution* Foam::faMesh::hasSolution() const
{
    return static_cast<const faSolution*>(this);
}


const Foam::faSchemes& Foam::faMesh::schemes() const
{
    return static_cast<const faSchemes&>(*this);
}


Foam::faSchemes& Foam::faMesh::schemes()
{
    return static_cast<faSchemes&>(*this);
}


const Foam::faSolution& Foam::faMesh::solution() const
{
    return static_cast<const faSolution&>(*this);
}


Foam::faSolution& Foam::faMesh::solution()
{
    return static_cast<faSolution&>(*this);
}


const Foam::polyMesh& Foam::faMesh::mesh() const
{
    return refCast<const polyMesh>(faMeshRegistry::parent().parent());
}


Foam::fileName Foam::faMesh::meshDir() const
{
    return dbDir()/faMesh::meshSubDir;
}


const Foam::Time& Foam::faMesh::time() const
{
    return faMeshRegistry::time();
}


const Foam::fileName& Foam::faMesh::pointsInstance() const
{
    return mesh().pointsInstance();
}


const Foam::fileName& Foam::faMesh::facesInstance() const
{
    return mesh().facesInstance();
}


const Foam::word& Foam::faMesh::regionName() const
{
    return polyMesh::regionName(objectRegistry::name());
}


Foam::labelList Foam::faMesh::faceCells() const
{
    const labelList& faceOwner = this->mesh().faceOwner();

    labelList list(faceLabels_);

    for (label& val : list)
    {
        // Transcribe from faceId to cellId (owner)
        val = faceOwner[val];
    }

    return list;
}


void Foam::faMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = thisDb().time().path()/instanceDir/meshDir();

    Foam::rm(meshFilesPath/"faceLabels");
    Foam::rm(meshFilesPath/"faBoundary");
}


void Foam::faMesh::removeFiles() const
{
    removeFiles(thisDb().instance());
}


const Foam::labelList& Foam::faMesh::patchStarts() const
{
    if (!patchStartsPtr_)
    {
        calcPatchStarts();
    }

    return *patchStartsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::Le() const
{
    if (!LePtr_)
    {
        calcLe();
    }

    return *LePtr_;
}


const Foam::edgeScalarField& Foam::faMesh::magLe() const
{
    if (!magLePtr_)
    {
        calcMagLe();
    }

    return *magLePtr_;
}


const Foam::areaVectorField& Foam::faMesh::areaCentres() const
{
    if (!faceCentresPtr_)
    {
        calcFaceCentres();
    }

    return *faceCentresPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeCentres() const
{
    if (!edgeCentresPtr_)
    {
        calcEdgeCentres();
    }

    return *edgeCentresPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S() const
{
    if (!SPtr_)
    {
        calcS();
    }

    return *SPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S0() const
{
    if (!S0Ptr_)
    {
        FatalErrorInFunction
            << "S0 is not available"
            << abort(FatalError);
    }

    return *S0Ptr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S00() const
{
    if (!S00Ptr_)
    {
        S00Ptr_ = std::make_unique<DimensionedField<scalar, areaMesh>>
        (
            IOobject
            (
                "S00",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            S0()
        );

        S0Ptr_->writeOpt(IOobject::AUTO_WRITE);
    }

    return *S00Ptr_;
}


const Foam::areaVectorField& Foam::faMesh::faceAreaNormals() const
{
    if (!faceAreaNormalsPtr_)
    {
        calcFaceAreaNormals();
    }

    return *faceAreaNormalsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeAreaNormals() const
{
    if (!edgeAreaNormalsPtr_)
    {
        calcEdgeAreaNormals();
    }

    return *edgeAreaNormalsPtr_;
}


const Foam::vectorField& Foam::faMesh::pointAreaNormals() const
{
    if (!pointAreaNormalsPtr_)
    {
        pointAreaNormalsPtr_ = std::make_unique<vectorField>(nPoints());

        calcPointAreaNormals(*pointAreaNormalsPtr_);

        if (quadricsFit_ > 0)
        {
            calcPointAreaNormalsByQuadricsFit(*pointAreaNormalsPtr_);
        }
    }

    return *pointAreaNormalsPtr_;
}


const Foam::areaScalarField& Foam::faMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        calcFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


const Foam::FieldField<Foam::Field, Foam::tensor>&
Foam::faMesh::edgeTransformTensors() const
{
    if (!edgeTransformTensorsPtr_)
    {
        calcEdgeTransformTensors();
    }

    return *edgeTransformTensorsPtr_;
}


bool Foam::faMesh::hasGlobalData() const noexcept
{
    return bool(globalMeshDataPtr_);
}


const Foam::faGlobalMeshData& Foam::faMesh::globalData() const
{
    if (!globalMeshDataPtr_)
    {
        globalMeshDataPtr_.reset(new faGlobalMeshData(*this));
    }

    return *globalMeshDataPtr_;
}


const Foam::lduAddressing& Foam::faMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        calcLduAddressing();
    }

    return *lduPtr_;
}


bool Foam::faMesh::movePoints()
{
    // Grab point motion from polyMesh
    const vectorField& newPoints = mesh().points();

    // Grab old time areas if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        if (S00Ptr_ && S0Ptr_)
        {
            DebugInfo<< "Copy old-old S" << endl;
            *S00Ptr_ = *S0Ptr_;
        }

        if (S0Ptr_)
        {
            DebugInfo<< "Copy old S" << endl;
            *S0Ptr_ = S();
        }
        else
        {
            DebugInfo<< "Creating old cell volumes." << endl;

            S0Ptr_ = std::make_unique<DimensionedField<scalar, areaMesh>>
            (
                IOobject
                (
                    "S0",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                S()
            );
        }

        curTimeIndex_ = time().timeIndex();
    }

    clearGeomNotAreas();

    if (patchPtr_)
    {
        patchPtr_->movePoints(newPoints);
    }

    // Move boundary points
    boundary_.movePoints(newPoints);

    // Move interpolation
    edgeInterpolation::movePoints();

    // Note: Fluxes were dummy?

    syncGeom();

    return true;
}


bool Foam::faMesh::correctPatchPointNormals(const label patchID) const
{
    if
    (
        bool(correctPatchPointNormalsPtr_)
     && patchID >= 0 && patchID < boundary().size()
    )
    {
        return (*correctPatchPointNormalsPtr_)[patchID];
    }

    return false;
}


Foam::boolList& Foam::faMesh::correctPatchPointNormals() const
{
    if (!correctPatchPointNormalsPtr_)
    {
        correctPatchPointNormalsPtr_ =
            std::make_unique<boolList>(boundary().size(), false);
    }

    return *correctPatchPointNormalsPtr_;
}


bool Foam::faMesh::write(const bool writeOnProc) const
{
    faceLabels_.write();
    boundary_.write();

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::faMesh::operator!=(const faMesh& m) const
{
    return &m != this;
}


bool Foam::faMesh::operator==(const faMesh& m) const
{
    return &m == this;
}


// ************************************************************************* //

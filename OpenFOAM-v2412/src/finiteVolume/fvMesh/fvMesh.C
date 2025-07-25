/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "fvMeshLduAddressing.H"
#include "mapPolyMesh.H"
#include "MapFvFields.H"
#include "fvMeshMapper.H"
#include "mapClouds.H"
#include "MeshObject.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::clearGeomNotOldVol()
{
    meshObject::clearUpto
    <
        fvMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    meshObject::clearUpto
    <
        lduMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    deleteDemandDrivenData(VPtr_);
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfPtr_);
}


void Foam::fvMesh::updateGeomNotOldVol()
{
    const bool haveV = (VPtr_ != nullptr);
    const bool haveSf = (SfPtr_ != nullptr);
    const bool haveMagSf = (magSfPtr_ != nullptr);
    const bool haveCP = (CPtr_ != nullptr);
    const bool haveCf = (CfPtr_ != nullptr);

    clearGeomNotOldVol();

    // Now recreate the fields
    if (haveV)
    {
        (void)V();
    }

    if (haveSf)
    {
        (void)Sf();
    }

    if (haveMagSf)
    {
        (void)magSf();
    }

    if (haveCP)
    {
        (void)C();
    }

    if (haveCf)
    {
        (void)Cf();
    }
}


void Foam::fvMesh::clearGeom()
{
    clearGeomNotOldVol();

    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);

    // Mesh motion flux cannot be deleted here because the old-time flux
    // needs to be saved.
}


void Foam::fvMesh::clearAddressing(const bool isMeshUpdate)
{
    DebugInFunction << "isMeshUpdate: " << isMeshUpdate << endl;

    if (isMeshUpdate)
    {
        // Part of a mesh update. Keep meshObjects that have an updateMesh
        // callback
        meshObject::clearUpto
        <
            fvMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
        meshObject::clearUpto
        <
            lduMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
    }
    else
    {
        meshObject::clear<fvMesh, TopologicalMeshObject>(*this);
        meshObject::clear<lduMesh, TopologicalMeshObject>(*this);
    }
    deleteDemandDrivenData(lduPtr_);
}


void Foam::fvMesh::storeOldVol(const scalarField& V)
{
    if (curTimeIndex_ < time().timeIndex())
    {
        DebugInFunction
            << " Storing old time volumes since from time " << curTimeIndex_
            << " and time now " << time().timeIndex()
            << " V:" << V.size() << endl;

        if (V00Ptr_ && V0Ptr_)
        {
            // Copy V0 into V00 storage
            *V00Ptr_ = *V0Ptr_;
        }

        if (!V0Ptr_)
        {
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "V0",
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                *this,
                dimVolume
            );
            // Note: V0 now sized with current mesh, not with (potentially
            //       different size) V.
            V0Ptr_->scalarField::resize_nocopy(V.size());
        }

        // Copy V into V0 storage
        V0Ptr_->scalarField::operator=(V);

        curTimeIndex_ = time().timeIndex();

        if (debug)
        {
            InfoInFunction
                << " Stored old time volumes V0:" << V0Ptr_->size()
                << endl;

            if (V00Ptr_)
            {
                InfoInFunction
                    << " Stored oldold time volumes V00:" << V00Ptr_->size()
                    << endl;
            }
        }
    }
}


void Foam::fvMesh::clearOutLocal(const bool isMeshUpdate)
{
    clearGeom();
    surfaceInterpolation::clearOut();

    clearAddressing(isMeshUpdate);

    // Clear mesh motion flux
    phiPtr_.reset(nullptr);
}


void Foam::fvMesh::clearOut(const bool isMeshUpdate)
{
    clearOutLocal(isMeshUpdate);

    polyMesh::clearOut(isMeshUpdate);
}


void Foam::fvMesh::clearMeshPhi()
{
    // Clear mesh motion flux
    phiPtr_.reset(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMesh::fvMesh(const IOobject& io, const bool doInit)
:
    polyMesh(io, doInit),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    surfaceInterpolation(*this),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    DebugInFunction << "Constructing fvMesh from IOobject" << endl;

    if (doInit)
    {
        fvMesh::init(false);    // do not initialise lower levels
    }
}


bool Foam::fvMesh::init(const bool doInit)
{
    if (doInit)
    {
        // Construct basic geometry calculation engine. Note: do before
        // doing anything with primitiveMesh::cellCentres etc.
        (void)geometry();

        // Initialise my data
        polyMesh::init(doInit);
    }

    // Read some optional fields
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // For redistributePar + fileHandler can have master processor
    // find the file but not the subprocs.

    IOobject rio
    (
        "none",
        this->time().timeName(),
        *this,
        IOobject::LAZY_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );


    // Read old cell volumes (if present) and set the storage of V00
    rio.resetHeader("V0");
    if (returnReduceOr(rio.typeHeaderOk<regIOobject>(false)))
    {
        DebugInFunction
            << "Detected V0: " << rio.objectRelPath() << nl;

        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            rio,
            *this,
            dimensionedScalar(dimVol, Foam::zero{})
        );

        // Set the moving flag early in case demand-driven geometry
        // construction checks for it
        moving(true);

        (void)V00();
    }


    // Read mesh fluxes (if present) and set the mesh to be moving
    rio.resetHeader("meshPhi");
    if (returnReduceOr(rio.typeHeaderOk<regIOobject>(false)))
    {
        DebugInFunction
            << "Detected meshPhi: " << rio.objectRelPath() << nl;

        // Clear mesh motion flux
        phiPtr_.reset(nullptr);

        phiPtr_ = std::make_unique<surfaceScalarField>
        (
            rio,
            *this,
            dimensionedScalar(dimVol/dimTime, Foam::zero{})
        );

        // Set the moving flag early in case demand-driven geometry
        // construction checks for it
        moving(true);

        // The mesh is now considered moving so the old-time cell volumes
        // will be required for the time derivatives so if they haven't been
        // read initialise to the current cell volumes
        if (!V0Ptr_)
        {
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "V0",
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                V()
            );
        }
    }

    // Assume something changed
    return true;
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    surfaceInterpolation(*this),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    DebugInFunction << "Constructing fvMesh from components" << endl;
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(cells),
        syncPar
    ),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    surfaceInterpolation(*this),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    DebugInFunction << "Constructing fvMesh from components" << endl;
}


Foam::fvMesh::fvMesh(const IOobject& io, const Foam::zero, const bool syncPar)
:
    fvMesh(io, pointField(), faceList(), labelList(), labelList(), syncPar)
{}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const fvMesh& baseMesh,
    const Foam::zero,
    const bool syncPar
)
:
    fvMesh
    (
        io,
        baseMesh,
        pointField(),
        faceList(),
        labelList(),  // owner
        labelList(),  // neighbour
        syncPar
    )
{}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const fvMesh& baseMesh,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour),
        syncPar
    ),
    fvSchemes
    (
        static_cast<const objectRegistry&>(*this),
        io.readOpt(),
        static_cast<const dictionary*>(baseMesh.hasSchemes())
    ),
    surfaceInterpolation(*this),
    fvSolution
    (
        static_cast<const objectRegistry&>(*this),
        io.readOpt(),
        static_cast<const dictionary*>(baseMesh.hasSolution())
    ),
    boundary_(*this),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    DebugInFunction << "Constructing fvMesh as copy and primitives" << endl;

    // Reset mesh data
    data().reset(baseMesh.data());
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const fvMesh& baseMesh,
    pointField&& points,
    faceList&& faces,
    cellList&& cells,
    const bool syncPar
)
:
    polyMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(cells),
        syncPar
    ),
    fvSchemes
    (
        static_cast<const objectRegistry&>(*this),
        io.readOpt(),
        static_cast<const dictionary*>(baseMesh.hasSchemes())
    ),
    surfaceInterpolation(*this),
    fvSolution
    (
        static_cast<const objectRegistry&>(*this),
        io.readOpt(),
        static_cast<const dictionary*>(baseMesh.hasSolution())
    ),
    boundary_(*this),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfPtr_(nullptr),
    magSfPtr_(nullptr),
    CPtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    DebugInFunction << "Constructing fvMesh as copy and primitives" << endl;

    // Reset mesh data
    data().reset(baseMesh.data());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvSchemes* Foam::fvMesh::hasSchemes() const
{
    return static_cast<const fvSchemes*>(this);
}


const Foam::fvSolution* Foam::fvMesh::hasSolution() const
{
    return static_cast<const fvSolution*>(this);
}


const Foam::fvSchemes& Foam::fvMesh::schemes() const
{
    return static_cast<const fvSchemes&>(*this);
}


Foam::fvSchemes& Foam::fvMesh::schemes()
{
    return static_cast<fvSchemes&>(*this);
}


const Foam::fvSolution& Foam::fvMesh::solution() const
{
    return static_cast<const fvSolution&>(*this);
}


Foam::fvSolution& Foam::fvMesh::solution()
{
    return static_cast<fvSolution&>(*this);
}


Foam::SolverPerformance<Foam::scalar> Foam::fvMesh::solve
(
    fvMatrix<scalar>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::vector> Foam::fvMesh::solve
(
    fvMatrix<vector>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::sphericalTensor> Foam::fvMesh::solve
(
    fvMatrix<sphericalTensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::symmTensor> Foam::fvMesh::solve
(
    fvMatrix<symmTensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::tensor> Foam::fvMesh::solve
(
    fvMatrix<tensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


void Foam::fvMesh::addFvPatches
(
    polyPatchList& plist,
    const bool validBoundary
)
{
    if (boundary().size())
    {
        FatalErrorInFunction
            << " boundary already exists"
            << abort(FatalError);
    }

    addPatches(plist, validBoundary);
    boundary_.addPatches(boundaryMesh());
}


void Foam::fvMesh::addFvPatches
(
    const List<polyPatch*>& p,
    const bool validBoundary
)
{
    // Acquire ownership of the pointers
    polyPatchList plist(const_cast<List<polyPatch*>&>(p));

    addFvPatches(plist, validBoundary);
}


void Foam::fvMesh::removeFvBoundary()
{
    DebugInFunction << "Removing boundary patches." << endl;

    // Remove fvBoundaryMesh data first.
    boundary_.clear();
    polyMesh::removeBoundary();

    clearOut();
}


Foam::polyMesh::readUpdateState Foam::fvMesh::readUpdate()
{
    DebugInFunction << "Updating fvMesh";

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        DebugInfo << "Boundary and topological update" << endl;

        boundary_.readUpdate(boundaryMesh());

        clearOut();

    }
    else if (state == polyMesh::TOPO_CHANGE)
    {
        DebugInfo << "Topological update" << endl;

        // fvMesh::clearOut() but without the polyMesh::clearOut
        clearOutLocal();
    }
    else if (state == polyMesh::POINTS_MOVED)
    {
        DebugInfo << "Point motion update" << endl;

        clearGeom();
    }
    else
    {
        DebugInfo << "No update" << endl;
    }

    return state;
}


const Foam::lduAddressing& Foam::fvMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        DebugInFunction
            << "Calculating fvMeshLduAddressing from nFaces:"
            << nFaces() << endl;

        lduPtr_ = new fvMeshLduAddressing(*this);

        return *lduPtr_;
    }

    return *lduPtr_;
}


Foam::lduInterfacePtrsList Foam::fvMesh::interfaces() const
{
    return boundary().interfaces();
}


void Foam::fvMesh::mapFields(const mapPolyMesh& meshMap)
{
    DebugInFunction
        << " nOldCells:" << meshMap.nOldCells()
        << " nCells:" << nCells()
        << " nOldFaces:" << meshMap.nOldFaces()
        << " nFaces:" << nFaces()
        << endl;

    // We require geometric properties valid for the old mesh
    if
    (
        meshMap.cellMap().size() != nCells()
     || meshMap.faceMap().size() != nFaces()
    )
    {
        FatalErrorInFunction
            << "mapPolyMesh does not correspond to the old mesh."
            << " nCells:" << nCells()
            << " cellMap:" << meshMap.cellMap().size()
            << " nOldCells:" << meshMap.nOldCells()
            << " nFaces:" << nFaces()
            << " faceMap:" << meshMap.faceMap().size()
            << " nOldFaces:" << meshMap.nOldFaces()
            << exit(FatalError);
    }

    // Create a mapper
    const fvMeshMapper mapper(*this, meshMap);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<vector, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<sphericalTensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<symmTensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<tensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<vector, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields
    <
        sphericalTensor, fvsPatchField, fvMeshMapper, surfaceMesh
    >
    (mapper);
    MapGeometricFields<symmTensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<tensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);

    // Map all the dimensionedFields in the objectRegistry
    MapDimensionedFields<scalar, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<vector, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<sphericalTensor, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<symmTensor, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<tensor, fvMeshMapper, volMesh>(mapper);

    // Map all the clouds in the objectRegistry
    mapClouds(*this, meshMap);


    const labelList& cellMap = meshMap.cellMap();

    // Map the old volume. Just map to new cell labels.
    if (V0Ptr_)
    {
        scalarField& V0 = *V0Ptr_;

        scalarField savedV0(V0);
        V0.resize_nocopy(nCells());

        forAll(V0, i)
        {
            if (cellMap[i] > -1)
            {
                V0[i] = savedV0[cellMap[i]];
            }
            else
            {
                V0[i] = 0.0;
            }
        }

        // Inject volume of merged cells
        label nMerged = 0;
        forAll(meshMap.reverseCellMap(), oldCelli)
        {
            label index = meshMap.reverseCellMap()[oldCelli];

            if (index < -1)
            {
                label celli = -index-2;

                V0[celli] += savedV0[oldCelli];

                nMerged++;
            }
        }

        DebugInfo
            << "Mapping old time volume V0. Merged "
            << nMerged << " out of " << nCells() << " cells" << endl;
    }


    // Map the old-old volume. Just map to new cell labels.
    if (V00Ptr_)
    {
        scalarField& V00 = *V00Ptr_;

        scalarField savedV00(V00);
        V00.resize_nocopy(nCells());

        forAll(V00, i)
        {
            if (cellMap[i] > -1)
            {
                V00[i] = savedV00[cellMap[i]];
            }
            else
            {
                V00[i] = 0.0;
            }
        }

        // Inject volume of merged cells
        label nMerged = 0;
        forAll(meshMap.reverseCellMap(), oldCelli)
        {
            label index = meshMap.reverseCellMap()[oldCelli];

            if (index < -1)
            {
                label celli = -index-2;

                V00[celli] += savedV00[oldCelli];
                nMerged++;
            }
        }

        DebugInfo
            << "Mapping old time volume V00. Merged "
            << nMerged << " out of " << nCells() << " cells" << endl;
    }
}


void Foam::fvMesh::movePoints(const pointField& p)
{
    DebugInFunction << endl;

    // Grab old time volumes if the time has been incremented
    // This will update V0, V00
    if (curTimeIndex_ < time().timeIndex())
    {
        storeOldVol(V());
    }


    // Move the polyMesh and initialise the mesh motion fluxes field
    // Note: mesh flux updated by the fvGeometryScheme

    if (!phiPtr_)
    {
        DebugInFunction<< "Creating initial meshPhi field" << endl;

        // Create mesh motion flux
        phiPtr_ = std::make_unique<surfaceScalarField>
        (
            IOobject
            (
                "meshPhi",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this,
            dimensionedScalar(dimVolume/dimTime, Foam::zero{})
        );
    }
    else
    {
        // Grab old time mesh motion fluxes if the time has been incremented
        if (phiPtr_->timeIndex() != time().timeIndex())
        {
            DebugInFunction<< "Accessing old-time meshPhi field" << endl;
            phiPtr_->oldTime();
        }
    }

    polyMesh::movePoints(p);

    // Update or delete the local geometric properties as early as possible so
    // they can be used if necessary. These get recreated here instead of
    // demand driven since they might do parallel transfers which can conflict
    // with when they're actually being used.
    // Note that between above "polyMesh::movePoints(p)" and here nothing
    // should use the local geometric properties.
    updateGeomNotOldVol();

    // Update other local data
    boundary_.movePoints();

    // Clear weights, deltaCoeffs, nonOrthoDeltaCoeffs, nonOrthCorrectionVectors
    surfaceInterpolation::clearOut();

    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);
}


void Foam::fvMesh::updateGeom()
{
    DebugInFunction << endl;

    // Let surfaceInterpolation handle geometry calculation. Note: this does
    // lower levels updateGeom
    surfaceInterpolation::updateGeom();
}


void Foam::fvMesh::updateMesh(const mapPolyMesh& mpm)
{
    DebugInFunction << endl;

    // Update polyMesh. This needs to keep volume existent!
    polyMesh::updateMesh(mpm);

    // Our slice of the addressing is no longer valid
    deleteDemandDrivenData(lduPtr_);

    if (VPtr_)
    {
        // Grab old time volumes if the time has been incremented
        // This will update V0, V00
        storeOldVol(mpm.oldCellVolumes());

        // Few checks
        if (VPtr_ && (VPtr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V:" << VPtr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
        if (V0Ptr_ && (V0Ptr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V0:" << V0Ptr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
        if (V00Ptr_ && (V00Ptr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V0:" << V00Ptr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
    }


    // Clear mesh motion flux (note: could instead save & map like volumes)
    if (phiPtr_)
    {
        // Mesh moving and topology change. Recreate meshPhi
        phiPtr_.reset(nullptr);

        // Create mesh motion flux
        phiPtr_ = std::make_unique<surfaceScalarField>
        (
            IOobject
            (
                "meshPhi",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *this,
            dimensionedScalar(dimVolume/dimTime, Foam::zero{})
        );
    }

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Map all fields
    mapFields(mpm);

    // Clear the current volume and other geometry factors
    surfaceInterpolation::updateMesh(mpm);

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObject::updateMesh<fvMesh>(*this, mpm);
    meshObject::updateMesh<lduMesh>(*this, mpm);
}


bool Foam::fvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok = true;
    if (phiPtr_)
    {
        ok = phiPtr_->write(writeOnProc);
        // NOTE: The old old time mesh phi might be necessary for certain
        // solver smooth restart using second order time schemes.
        //ok = phiPtr_->oldTime().write();
    }
    if (V0Ptr_ && V0Ptr_->writeOpt() == IOobject::AUTO_WRITE)
    {
        // For second order restarts we need to write V0
        ok = V0Ptr_->write(writeOnProc);
    }

    return ok && polyMesh::writeObject(streamOpt, writeOnProc);
}


bool Foam::fvMesh::write(const bool writeOnProc) const
{
    return polyMesh::write(writeOnProc);
}


template<>
typename Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::fvMesh::validComponents<Foam::sphericalTensor>() const
{
    return Foam::pTraits<Foam::sphericalTensor>::labelType(1);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::fvMesh::operator!=(const fvMesh& rhs) const
{
    return &rhs != this;
}


bool Foam::fvMesh::operator==(const fvMesh& rhs) const
{
    return &rhs == this;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "MapFieldConstraint.H"
#include "fvMatrices.H"
#include "meshToMesh.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static inline tmp<volScalarField> createField
(
    const fvMesh& mesh,
    const scalar val
)
{
    return volScalarField::New
    (
        polyMesh::defaultRegion,
        IOobject::NO_REGISTER,
        mesh,
        val,
        dimless
    );
}

}  // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::MapFieldConstraint<Type>::setSourceMesh
(
    refPtr<fvMesh>& meshRef,
    const autoPtr<Time>& runTimePtr
)
{
    const Time& runTime = runTimePtr();
    const word meshName(polyMesh::defaultRegion);

    // Fetch mesh from Time database
    meshRef.cref
    (
        runTime.cfindObject<fvMesh>(meshName)
    );

    if (!meshRef)
    {
        // Fallback: load mesh from disk and cache it
        meshRef.reset
        (
            new fvMesh
            (
                IOobject
                (
                    meshName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                )
            )
        );
    }
}


template<class Type>
void Foam::fv::MapFieldConstraint<Type>::createInterpolation
(
    const fvMesh& srcMesh,
    const fvMesh& tgtMesh
)
{
    if (consistent_)
    {
        interpPtr_.reset
        (
            new meshToMesh
            (
                srcMesh,
                tgtMesh,
                mapMethodName_,
                patchMapMethodName_
            )
        );
    }
    else
    {
        interpPtr_.reset
        (
            new meshToMesh
            (
                srcMesh,
                tgtMesh,
                mapMethodName_,
                patchMapMethodName_,
                patchMap_,
                cuttingPatches_
            )
        );
    }
}


template<class Type>
template<class VolFieldType>
VolFieldType& Foam::fv::MapFieldConstraint<Type>::getOrReadField
(
    const fvMesh& thisMesh,
    const word& fieldName
) const
{
    auto* ptr = thisMesh.getObjectPtr<VolFieldType>(fieldName);

    if (!ptr)
    {
        ptr = new VolFieldType
        (
            IOobject
            (
                fieldName,
                thisMesh.time().timeName(),
                thisMesh.thisDb(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            ),
            thisMesh
        );
        regIOobject::store(ptr);
    }

    return *ptr;
}


template<class Type>
Foam::labelList Foam::fv::MapFieldConstraint<Type>::tgtCellIDs() const
{
    const fvMesh& srcMesh = srcMeshPtr_();
    const fvMesh& tgtMesh = mesh_;

    // Create mask fields
    const volScalarField srcFld(createField(srcMesh, 1));
    volScalarField tgtFld(createField(tgtMesh, 0));

    // Map the mask field of 1s onto the mask field of 0s
    interpPtr_->mapSrcToTgt(srcFld, plusEqOp<scalar>(), tgtFld);

    // Identify and collect cell indices whose values were changed from 0 to 1
    DynamicList<label> cells;
    forAll(tgtFld, i)
    {
        if (tgtFld[i] != 0)
        {
            cells.append(i);
        }
    }

    return cells;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::MapFieldConstraint<Type>::transform::transform()
:
    positionPtr_(),
    directionPtr_(),
    points_(),
    origin_(),
    normal_(),
    active_(false)
{}


template<class Type>
Foam::fv::MapFieldConstraint<Type>::MapFieldConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    transform_(),
    srcTimePtr_(),
    srcMeshPtr_(),
    interpPtr_(),
    patchMap_(),
    cells_(),
    cuttingPatches_(),
    mapMethodName_(),
    patchMapMethodName_(),
    consistent_(false)
{
    read(dict);

    setSourceMesh(srcMeshPtr_, srcTimePtr_);

    const fvMesh& srcMesh = srcMeshPtr_();
    const fvMesh& tgtMesh = mesh_;
    createInterpolation(srcMesh, tgtMesh);

    cells_ = tgtCellIDs();

    if (returnReduceAnd(cells_.empty()))
    {
        WarningInFunction
            << "No cells selected!" << endl;
    }

    transform_.initialize(srcMesh, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::MapFieldConstraint<Type>::transform::initialize
(
    const fvMesh& srcMesh,
    const dictionary& dict
)
{
    const dictionary* subDictPtr = dict.findDict("transform");

    if (!subDictPtr)
    {
        return false;
    }

    positionPtr_.reset
    (
        Function1<point>::NewIfPresent
        (
            "position",
            *subDictPtr,
            word::null,
            &srcMesh
        )
    );

    directionPtr_.reset
    (
        Function1<point>::NewIfPresent
        (
            "direction",
            *subDictPtr,
            word::null,
            &srcMesh
        )
    );

    if (positionPtr_)
    {
        subDictPtr->readIfPresent("origin", origin_);
    }

    if (directionPtr_)
    {
        subDictPtr->readIfPresent("normal", normal_);
        normal_.normalise();
    }

    points_ = srcMesh.points();

    active_ = true;

    return true;
}


template<class Type>
void Foam::fv::MapFieldConstraint<Type>::transform::translate
(
    refPtr<fvMesh>& srcMeshPtr,
    const scalar t
)
{
    if (!positionPtr_)
    {
        return;
    }

    const pointField translate
    (
        points_ + (positionPtr_->value(t) - origin_)
    );

    fvMesh& srcMesh = srcMeshPtr.ref();
    srcMesh.movePoints(translate);
}


template<class Type>
void Foam::fv::MapFieldConstraint<Type>::transform::rotate
(
    refPtr<fvMesh>& srcMeshPtr,
    const scalar t
)
{
    if (!directionPtr_)
    {
        return;
    }

    const vector dir(normalised(directionPtr_->value(t)));

    const tensor rot(rotationTensor(normal_, dir));

    pointField rotate(points_);

    Foam::transform(rotate, rot, rotate);

    fvMesh& srcMesh = srcMeshPtr.ref();
    srcMesh.movePoints(rotate);
}


template<class Type>
bool Foam::fv::MapFieldConstraint<Type>::read(const dictionary& dict)
{
    if (!fv::option::read(dict))
    {
        return false;
    }

    fieldNames_.resize(1, coeffs_.getWord("field"));

    fv::option::resetApplied();

    // Load the time database for the source mesh once per simulation
    if (!srcTimePtr_)
    {
        fileName srcMesh(coeffs_.get<fileName>("srcMesh").expand());
        srcMesh.clean();

        srcTimePtr_.reset(Time::New(srcMesh));

        // Set time-step of source database to an arbitrary yet safe value
        srcTimePtr_().setDeltaT(1.0);
    }

    coeffs_.readEntry("mapMethod", mapMethodName_);
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName_))
    {
        FatalIOErrorInFunction(coeffs_)
            << type() << " " << name() << ": unknown map method "
            << mapMethodName_ << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_
            << exit(FatalIOError);
    }

    coeffs_.readIfPresent("consistent", consistent_);
    coeffs_.readIfPresent("patchMap", patchMap_);
    coeffs_.readIfPresent("cuttingPatches", cuttingPatches_);

    if (!coeffs_.readIfPresent("patchMapMethod", patchMapMethodName_))
    {
        meshToMesh::interpolationMethod mapMethod
        (
            meshToMesh::interpolationMethodNames_[mapMethodName_]
        );
        patchMapMethodName_ = meshToMesh::interpolationMethodAMI(mapMethod);
    }

    return true;
}


template<class Type>
void Foam::fv::MapFieldConstraint<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label
)
{
    DebugInfo
        << "MapFieldConstraint<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    // Translate and/or rotate source mesh if requested
    if (transform_.isActive())
    {
        // Use time from mesh_ since source mesh does not advance in time
        const scalar t = mesh_.time().timeOutputValue();
        transform_.translate(srcMeshPtr_, t);
        transform_.rotate(srcMeshPtr_, t);
    }

    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const word& fldName = fieldNames_[0];

    const fvMesh& srcMesh = srcMeshPtr_();
    const fvMesh& tgtMesh = mesh_;

    // Fetch source and target fields
    const VolFieldType& srcFld = getOrReadField<VolFieldType>(srcMesh, fldName);
    VolFieldType& tgtFld = tgtMesh.lookupObjectRef<VolFieldType>(fldName);

    // When mesh/src changes, reinitilize mesh-to-mesh members
    if (tgtMesh.changing() || transform_.isActive())
    {
        createInterpolation(srcMesh, tgtMesh);
        cells_ = tgtCellIDs();
    }

    // Map source-mesh field onto target-mesh field
    interpPtr_->mapSrcToTgt(srcFld, plusEqOp<Type>(), tgtFld);

    // Constrain mapped field in target mesh to avoid overwrite by solver
    eqn.setValues(cells_, UIndirectList<Type>(tgtFld, cells_));
}


// ************************************************************************* //

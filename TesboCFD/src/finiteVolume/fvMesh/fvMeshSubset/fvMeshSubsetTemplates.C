/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvMeshSubset.H"
#include "emptyFvsPatchField.H"
#include "emptyPointPatchField.H"
#include "emptyFvPatchFields.H"
#include "directFvPatchFieldMapper.H"
#include "directPointPatchFieldMapper.H"
#include "flipOp.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
>
Foam::fvMeshSubset::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const fvMesh& sMesh,
    const labelUList& patchMap,
    const labelUList& cellMap,
    const labelUList& faceMap,
    const bool allowUnmapped
)
{
    // 1. Create the complete field with dummy patch fields
    PtrList<fvPatchField<Type>> patchFields(patchMap.size());

    forAll(patchFields, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces. Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] == -1)
        {
            patchFields.set
            (
                patchi,
                new emptyFvPatchField<Type>
                (
                    sMesh.boundary()[patchi],
                    fvPatchField<Type>::Internal::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    fvPatchFieldBase::calculatedType(),
                    sMesh.boundary()[patchi],
                    fvPatchField<Type>::Internal::null()
                )
            );
        }
    }

    auto tresult = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "subset"+vf.name(),
            sMesh.time().timeName(),
            sMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sMesh,
        vf.dimensions(),
        Field<Type>(vf.primitiveField(), cellMap),
        patchFields
    );
    auto& result = tresult.ref();
    result.oriented() = vf.oriented();


    // 2. Change the fvPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(bf, patchi)
    {
        const label basePatchId = patchMap[patchi];

        if (basePatchId != -1)
        {
            // Construct addressing
            const fvPatch& subPatch = sMesh.boundary()[patchi];
            const fvPatch& basePatch = vf.mesh().boundary()[basePatchId];
            const label baseStart = basePatch.start();
            const label baseSize = basePatch.size();

            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                const label baseFacei = faceMap[subPatch.start()+i];

                if (baseFacei >= baseStart && baseFacei < baseStart+baseSize)
                {
                    directAddressing[i] = baseFacei-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Leave up to
                    // fvPatchField
                    directAddressing[i] = -1;
                }
            }


            directFvPatchFieldMapper mapper(directAddressing);

            // allowUnmapped : special mode for if we do not want to be
            // warned for unmapped faces (e.g. from fvMeshDistribute).

            const bool hasUnmapped = mapper.hasUnmapped();
            if (allowUnmapped)
            {
                mapper.hasUnmapped() = false;
            }

            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    vf.boundaryField()[basePatchId],
                    subPatch,
                    result.internalField(),
                    mapper
                )
            );

            if (allowUnmapped && hasUnmapped)
            {
                // Set unmapped values to zeroGradient. This is the default
                // action for unmapped fvPatchFields. Note that this bypasses
                // any special logic for handling unmapped fvPatchFields but
                // since this is only used inside fvMeshDistribute ...

                tmp<Field<Type>> tfld(bf[patchi].patchInternalField());
                const Field<Type>& fld = tfld();

                Field<Type> value(bf[patchi]);
                forAll(directAddressing, i)
                {
                    if (directAddressing[i] == -1)
                    {
                        value[i] = fld[i];
                    }
                }
                bf[patchi].fvPatchField<Type>::operator=(value);
            }
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>
>
Foam::fvMeshSubset::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const bool allowUnmapped
) const
{
    return interpolate
    (
        vf,
        subMesh(),
        patchMap(),
        cellMap(),
        faceMap(),
        allowUnmapped
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>
>
Foam::fvMeshSubset::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& vf,
    const fvMesh& sMesh,
    const labelUList& patchMap,
    const labelUList& cellMap,
    const labelUList& faceMap
)
{
    // 1. Create the complete field with dummy patch fields
    PtrList<fvsPatchField<Type>> patchFields(patchMap.size());

    forAll(patchFields, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces. Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] == -1)
        {
            patchFields.set
            (
                patchi,
                new emptyFvsPatchField<Type>
                (
                    sMesh.boundary()[patchi],
                    fvsPatchField<Type>::Internal::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    fvsPatchFieldBase::calculatedType(),
                    sMesh.boundary()[patchi],
                    fvsPatchField<Type>::Internal::null()
                )
            );
        }
    }

    // Create the complete field from the pieces
    auto tresult = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
    (
        IOobject
        (
            "subset"+vf.name(),
            sMesh.time().timeName(),
            sMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sMesh,
        vf.dimensions(),
        Field<Type>
        (
            vf.primitiveField(),
            SubList<label>(faceMap, sMesh.nInternalFaces())
        ),
        patchFields
    );
    auto& result = tresult.ref();
    result.oriented() = vf.oriented();


    // 2. Change the fvsPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(bf, patchi)
    {
        if (patchMap[patchi] != -1)
        {
            // Construct addressing
            const fvPatch& subPatch = sMesh.boundary()[patchi];
            const fvPatch& basePatch = vf.mesh().boundary()[patchMap[patchi]];
            const label baseStart = basePatch.start();
            const label baseSize = basePatch.size();

            labelList directAddressing(subPatch.size());

            forAll(directAddressing, i)
            {
                label baseFacei = faceMap[subPatch.start()+i];

                if (baseFacei >= baseStart && baseFacei < baseStart+baseSize)
                {
                    directAddressing[i] = baseFacei-baseStart;
                }
                else
                {
                    // Mapped from internal face. Do what? Leave up to
                    // patchField. This would require also to pass in
                    // original internal field so for now do as post-processing
                    directAddressing[i] = -1;
                }
            }

            bf.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchi]],
                    subPatch,
                    result(),
                    directFvPatchFieldMapper(directAddressing)
                )
            );


            // Post-process patch field for exposed faces

            fvsPatchField<Type>& pfld = bf[patchi];
            const labelUList& fc = bf[patchi].patch().faceCells();
            const labelList& own = vf.mesh().faceOwner();

            forAll(pfld, i)
            {
                label baseFacei = faceMap[subPatch.start()+i];
                if (baseFacei < vf.primitiveField().size())
                {
                    Type val = vf.internalField()[baseFacei];

                    if (cellMap[fc[i]] == own[baseFacei] || !vf.is_oriented())
                    {
                        pfld[i] = val;
                    }
                    else
                    {
                        pfld[i] = flipOp()(val);
                    }
                }
                else
                {
                    // Exposed face from other patch.
                    // Only possible in case of a coupled boundary
                    label patchi = vf.mesh().boundaryMesh().whichPatch
                    (
                        baseFacei
                    );
                    const fvPatch& otherPatch = vf.mesh().boundary()[patchi];
                    label patchFacei = otherPatch.patch().whichFace(baseFacei);
                    pfld[i] = vf.boundaryField()[patchi][patchFacei];
                }
            }
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>
>
Foam::fvMeshSubset::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sf,
    const bool allowUnmapped
) const
{
    return interpolate
    (
        sf,
        subMesh(),
        patchMap(),
        cellMap(),
        faceMap()
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>
>
Foam::fvMeshSubset::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& vf,
    const pointMesh& sMesh,
    const labelUList& patchMap,
    const labelUList& pointMap
)
{
    // 1. Create the complete field with dummy patch fields
    PtrList<pointPatchField<Type>> patchFields(patchMap.size());

    forAll(patchFields, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] == -1)
        {
            patchFields.set
            (
                patchi,
                new emptyPointPatchField<Type>
                (
                    sMesh.boundary()[patchi],
                    pointPatchField<Type>::Internal::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    pointPatchFieldBase::calculatedType(),
                    sMesh.boundary()[patchi],
                    pointPatchField<Type>::Internal::null()
                )
            );
        }
    }

    // Create the complete field from the pieces
    auto tresult = tmp<GeometricField<Type, pointPatchField, pointMesh>>::New
    (
        IOobject
        (
            "subset"+vf.name(),
            sMesh.time().timeName(),
            sMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sMesh,
        vf.dimensions(),
        Field<Type>(vf.primitiveField(), pointMap),
        patchFields
    );
    auto& result = tresult.ref();
    result.oriented() = vf.oriented();


    // 2. Change the pointPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(bf, patchi)
    {
        // Set the first one by hand as it corresponds to the
        // exposed internal faces.  Additional interpolation can be put here
        // as necessary.
        if (patchMap[patchi] != -1)
        {
            // Construct addressing
            const pointPatch& basePatch =
                vf.mesh().boundary()[patchMap[patchi]];

            const labelList& meshPoints = basePatch.meshPoints();

            // Make addressing from mesh to patch point
            Map<label> meshPointMap(invertToMap(meshPoints));


            // Find which subpatch points originate from which patch point
            const pointPatch& subPatch = sMesh.boundary()[patchi];
            const labelList& subMeshPoints = subPatch.meshPoints();

            // If mapped from outside patch leave handling up to patchField
            labelList directAddressing(subPatch.size(), -1);

            forAll(subMeshPoints, localI)
            {
                // Get mesh point on original mesh.
                label meshPointi = pointMap[subMeshPoints[localI]];

                const auto iter = meshPointMap.cfind(meshPointi);

                if (iter.good())
                {
                    directAddressing[localI] = iter.val();
                }
            }

            bf.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    vf.boundaryField()[patchMap[patchi]],
                    subPatch,
                    result.internalField(),
                    directPointPatchFieldMapper(directAddressing)
                )
            );
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>
>
Foam::fvMeshSubset::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& sf,
    const bool allowUnmapped
) const
{
    return interpolate
    (
        sf,
        pointMesh::New(subMesh()),     // subsetted point mesh
        pointPatchMap(),
        pointMap()
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::DimensionedField<Type, Foam::volMesh>
>
Foam::fvMeshSubset::interpolate
(
    const DimensionedField<Type, volMesh>& df,
    const fvMesh& sMesh,
    const labelUList& cellMap
)
{
    auto tresult = tmp<DimensionedField<Type, volMesh>>::New
    (
        IOobject
        (
            "subset"+df.name(),
            sMesh.time().timeName(),
            sMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sMesh,
        df.dimensions(),
        Field<Type>(df, cellMap)
    );
    tresult.ref().oriented() = df.oriented();

    return tresult;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::DimensionedField<Type, Foam::volMesh>
>
Foam::fvMeshSubset::interpolate
(
    const DimensionedField<Type, volMesh>& df,
    const bool allowUnmapped
) const
{
    return interpolate(df, subMesh(), cellMap());
}


// ************************************************************************* //

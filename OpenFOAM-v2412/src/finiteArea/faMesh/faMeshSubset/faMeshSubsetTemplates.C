/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "faMeshSubset.H"
#include "areaFaMesh.H"
#include "edgeFaMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "emptyFaPatchFields.H"
#include "directFaPatchFieldMapper.H"
#include "flipOp.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>
>
Foam::faMeshSubset::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const faMesh& sMesh,
    const bool allowUnmapped
)
{
    // 1. Create the complete field with dummy patch fields
    PtrList<faPatchField<Type>> patchFields(sMesh.boundary().size());

    forAll(patchFields, patchi)
    {
        patchFields.set
        (
            patchi,
            faPatchField<Type>::New
            (
                faPatchFieldBase::calculatedType(),
                sMesh.boundary()[patchi],
                faPatchField<Type>::Internal::null()
            )
        );
    }

    auto tresult = tmp<GeometricField<Type, faPatchField, areaMesh>>::New
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
        Field<Type>(),
        // Field<Type>(vf.primitiveField(), cellMap),
        patchFields
    );
    auto& result = tresult.ref();
    result.oriented() = vf.oriented();


    // 2. Change the faPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(bf, patchi)
    {
        // Construct addressing
        const faPatch& subPatch = sMesh.boundary()[patchi];

        labelList directAddressing;
        directFaPatchFieldMapper mapper(directAddressing);

        // allowUnmapped : special mode for if we do not want to be
        // warned for unmapped faces (e.g. from faMeshDistribute).
        const bool hasUnmapped = mapper.hasUnmapped();
        if (allowUnmapped)
        {
            mapper.hasUnmapped() = false;
        }

        bf.set
        (
            patchi,
            faPatchField<Type>::New
            (
                vf.boundaryField()[patchi],
                subPatch,
                result(),
                mapper
            )
        );

        if (allowUnmapped && hasUnmapped)
        {
            // Set unmapped values to zeroGradient. This is the default
            // action for unmapped faPatchFields. Note that this bypasses
            // any special logic for handling unmapped faPatchFields but
            // since this is only used inside faMeshDistribute ...

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
            bf[patchi].faPatchField<Type>::operator=(value);
        }
    }

    return tresult;
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>
>
Foam::faMeshSubset::interpolate
(
    const GeometricField<Type, faePatchField, edgeMesh>& vf,
    const faMesh& sMesh
)
{
    // 1. Create the complete field with dummy patch fields
    PtrList<faePatchField<Type>> patchFields(sMesh.boundary().size());

    forAll(patchFields, patchi)
    {
        patchFields.set
        (
            patchi,
            faePatchField<Type>::New
            (
                faePatchFieldBase::calculatedType(),
                sMesh.boundary()[patchi],
                faePatchField<Type>::Internal::null()
            )
        );
    }

    auto tresult = tmp<GeometricField<Type, faePatchField, edgeMesh>>::New
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
        Field<Type>(),
        // Field<Type>
        // (
        //     vf.primitiveField(),
        //     SubList<label>(edgeMap, sMesh.nInternalEdges())
        // ),
        patchFields
    );
    auto& result = tresult.ref();
    result.oriented() = vf.oriented();


    // 2. Change the faePatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(bf, patchi)
    {
        // Construct addressing
        const faPatch& subPatch = sMesh.boundary()[patchi];

        labelList directAddressing;
        directFaPatchFieldMapper mapper(directAddressing);

        bf.set
        (
            patchi,
            faePatchField<Type>::New
            (
                vf.boundaryField()[patchi],
                subPatch,
                result(),
                mapper
            )
        );
    }

    return tresult;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>
>
Foam::faMeshSubset::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const bool allowUnmapped
) const
{
    if (subMeshPtr_)
    {
        return interpolate(vf, *subMeshPtr_);
    }

    return vf;
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>
>
Foam::faMeshSubset::interpolate
(
    const GeometricField<Type, faePatchField, edgeMesh>& vf,
    const bool allowUnmapped
) const
{
    if (subMeshPtr_)
    {
        return interpolate(vf, *subMeshPtr_);
    }

    return vf;
}


// ************************************************************************* //

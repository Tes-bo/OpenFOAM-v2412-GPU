/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "fvFieldDecomposer.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"
#include "processorCyclicFvPatchField.H"
#include "processorCyclicFvsPatchField.H"
#include "emptyFvPatchFields.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldDecomposer::decomposeField
(
    const DimensionedField<Type, volMesh>& field
) const
{
    // Create the field for the processor

    return DimensionedField<Type, volMesh>::New
    (
        field.name(),
        IOobject::NO_REGISTER,
        procMesh_,
        field.dimensions(),
        // Internal field - mapped values
        Field<Type>(field.field(), cellAddressing_)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const bool allowUnknownPatchFields
) const
{
    // Create the field for the processor
    // - with dummy patch fields
    auto tresult = GeometricField<Type, fvPatchField, volMesh>::New
    (
        field.name(),
        IOobject::NO_REGISTER,
        procMesh_,
        field.dimensions(),
        // Internal field - mapped values
        Field<Type>(field.primitiveField(), cellAddressing_),
        fvPatchFieldBase::calculatedType()
    );
    auto& result = tresult.ref();
    result.oriented() = field.oriented();


    // 2. Change the fvPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(bf, patchi)
    {
        if (patchFieldDecomposerPtrs_.set(patchi))
        {
            bf.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    result.internalField(),
                    patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorCyclicFvPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    result.internalField(),
                    Field<Type>
                    (
                        field.primitiveField(),
                        processorVolPatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
        else if (isA<processorFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorFvPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    result.internalField(),
                    Field<Type>
                    (
                        field.primitiveField(),
                        processorVolPatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
        else if (allowUnknownPatchFields)
        {
            bf.set
            (
                patchi,
                new emptyFvPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    result.internalField()
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unknown type." << abort(FatalError);
        }
    }

    // Create the field for the processor
    return tresult;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    labelList mapAddr
    (
        labelList::subList(faceAddressing_, procMesh_.nInternalFaces())
    );
    forAll(mapAddr, i)
    {
        mapAddr[i] -= 1;
    }


    // Problem with addressing when a processor patch picks up both internal
    // faces and faces from cyclic boundaries. This is a bit of a hack, but
    // I cannot find a better solution without making the internal storage
    // mechanism for surfaceFields correspond to the one of faces in polyMesh
    // (i.e. using slices)
    Field<Type> allFaceField(field.mesh().nFaces());

    {
        SubList<Type>(allFaceField, field.primitiveField().size()) =
            field.primitiveField();

        forAll(field.boundaryField(), patchi)
        {
            const Field<Type>& pfld = field.boundaryField()[patchi];

            const label start = field.mesh().boundaryMesh()[patchi].start();

            SubList<Type>(allFaceField, pfld.size(), start) = pfld;
        }
    }


    // 1. Create the complete field with dummy patch fields

    auto tresult = GeometricField<Type, fvsPatchField, surfaceMesh>::New
    (
        field.name(),
        IOobject::NO_REGISTER,
        procMesh_,
        field.dimensions(),
        // Internal field - mapped values
        Field<Type>(field.primitiveField(), mapAddr),
        fvsPatchFieldBase::calculatedType()
    );
    auto& result = tresult.ref();
    result.oriented() = field.oriented();

    // 2. Change the fvsPatchFields to the correct type using a mapper
    //  constructor (with reference to the now correct internal field)

    auto& bf = result.boundaryFieldRef();

    forAll(boundaryAddressing_, patchi)
    {
        if (patchFieldDecomposerPtrs_.set(patchi))
        {
            bf.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchi]],
                    procMesh_.boundary()[patchi],
                    result.internalField(),
                    patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorCyclicFvsPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    result.internalField(),
                    Field<Type>
                    (
                        allFaceField,
                        processorSurfacePatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );

            if (result.is_oriented())
            {
                bf[patchi] *= faceSign_[patchi];
            }
        }
        else if (isA<processorFvPatch>(procMesh_.boundary()[patchi]))
        {
            bf.set
            (
                patchi,
                new processorFvsPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    result.internalField(),
                    Field<Type>
                    (
                        allFaceField,
                        processorSurfacePatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );

            if (result.is_oriented())
            {
                bf[patchi] *= faceSign_[patchi];
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unknown type." << abort(FatalError);
        }
    }

    // Create the field for the processor
    return tresult;
}


template<class GeoField>
void Foam::fvFieldDecomposer::decomposeFields
(
    const PtrList<GeoField>& fields
) const
{
    for (const auto& fld : fields)
    {
        decomposeField(fld)().write();
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "faFieldDecomposer.H"
#include "GeometricField.H"
#include "IOobjectList.H"
#include "processorFaPatchField.H"
#include "processorFaePatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>>
Foam::faFieldDecomposer::decomposeField
(
    const GeometricField<Type, faPatchField, areaMesh>& field
) const
{
    // Create and map the internal field values
    Field<Type> internalField(field.internalField(), faceAddressing_);

    // Create and map the patch field values
    PtrList<faPatchField<Type>> patchFields(boundaryAddressing_.size());

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];

        if (oldPatchi >= 0)
        {
            patchFields.set
            (
                patchi,
                faPatchField<Type>::New
                (
                    field.boundaryField()[oldPatchi],
                    procMesh_.boundary()[patchi],
                    faPatchField<Type>::Internal::null(),
                    patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                new processorFaPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    faPatchField<Type>::Internal::null(),
                    Field<Type>
                    (
                        field.internalField(),
                        processorAreaPatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
    }

    // Create the field for the processor
    return
        tmp<GeometricField<Type, faPatchField, areaMesh>>::New
        (
            IOobject
            (
                field.name(),
                procMesh_.thisDb().time().timeName(),
                procMesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            field.dimensions(),
            internalField,
            patchFields
        );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::faFieldDecomposer::decomposeField
(
    const GeometricField<Type, faePatchField, edgeMesh>& field
) const
{
    labelList mapAddr
    (
        labelList::subList
        (
            edgeAddressing_,
            procMesh_.nInternalEdges()
        )
    );
    forAll(mapAddr, i)
    {
        mapAddr[i] -= 1;
    }

    // Create and map the internal field values
    Field<Type> internalField
    (
        field.internalField(),
        mapAddr
    );

    // Problem with addressing when a processor patch picks up both internal
    // edges and edges from cyclic boundaries. This is a bit of a hack, but
    // I cannot find a better solution without making the internal storage
    // mechanism for edgeFields correspond to the one of edges in polyMesh
    // (i.e. using slices)
    Field<Type> allEdgeField(field.mesh().nEdges());

    forAll(field.internalField(), i)
    {
        allEdgeField[i] = field.internalField()[i];
    }

    forAll(field.boundaryField(), patchi)
    {
        const Field<Type>& p = field.boundaryField()[patchi];

        const label patchStart = field.mesh().boundary()[patchi].start();

        forAll(p, i)
        {
            allEdgeField[patchStart + i] = p[i];
        }
    }

    // Create and map the patch field values
    PtrList<faePatchField<Type>> patchFields(boundaryAddressing_.size());

    forAll(boundaryAddressing_, patchi)
    {
        const label oldPatchi = boundaryAddressing_[patchi];

        if (oldPatchi >= 0)
        {
            patchFields.set
            (
                patchi,
                faePatchField<Type>::New
                (
                    field.boundaryField()[oldPatchi],
                    procMesh_.boundary()[patchi],
                    faePatchField<Type>::Internal::null(),
                    patchFieldDecomposerPtrs_[patchi]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                new processorFaePatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    faePatchField<Type>::Internal::null(),
                    Field<Type>
                    (
                        allEdgeField,
                        processorEdgePatchFieldDecomposerPtrs_[patchi]
                    )
                )
            );
        }
    }

    // Create the field for the processor
    return
        tmp<GeometricField<Type, faePatchField, edgeMesh>>::New
        (
            IOobject
            (
                field.name(),
                procMesh_.thisDb().time().timeName(),
                procMesh_.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            field.dimensions(),
            internalField,
            patchFields
        );
}


template<class GeoField>
void Foam::faFieldDecomposer::decomposeFields
(
    const PtrList<GeoField>& fields
) const
{
    forAll(fields, fieldi)
    {
        decomposeField(fields[fieldi])().write();
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "pointFieldReconstructor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::pointFieldReconstructor::reconstructField
(
    const IOobject& fieldObject,
    const PtrList<GeometricField<Type, pointPatchField, pointMesh>>& procFields
) const
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    // Create the internalField
    Field<Type> internalField(mesh_.size());

    // Create the patch fields
    PtrList<pointPatchField<Type>> patchFields(mesh_.boundary().size());


    forAll(procMeshes_, proci)
    {
        const auto& procField = procFields[proci];

        // Get processor-to-global addressing for use in rmap
        const labelList& procToGlobalAddr = pointProcAddressing_[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            procToGlobalAddr
        );

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[proci], patchi)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                if (!patchFields.set(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        pointPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            pointPatchField<Type>::Internal::null(),
                            pointPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    patchPointAddressing_[proci][patchi]
                );
            }
        }
    }

    // Construct and write the field
    // setting the internalField and patchFields
    return tmp<fieldType>::New
    (
        IOobject
        (
            fieldObject.name(),
            mesh_.thisDb().time().timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        procFields[0].dimensions(),
        internalField,
        patchFields
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::pointFieldReconstructor::reconstructPointField
(
    const IOobject& fieldObject
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    // Read the field for all the processors
    PtrList<fieldType> procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new fieldType
            (
                IOobject
                (
                    fieldObject.name(),
                    procMeshes_[proci].thisDb().time().timeName(),
                    procMeshes_[proci].thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci]
            )
        );
    }

    return reconstructField
    (
        IOobject
        (
            fieldObject.name(),
            mesh_.thisDb().time().timeName(),
            mesh_.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        procFields
    );
}


template<class Type>
Foam::label Foam::pointFieldReconstructor::reconstructPointFields
(
    const UPtrList<const IOobject>& fieldObjects
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    label nFields = 0;

    for (const IOobject& io : fieldObjects)
    {
        if (io.isHeaderClass<fieldType>())
        {
            if (verbose_)
            {
                if (!nFields)
                {
                    Info<< "    Reconstructing "
                        << fieldType::typeName << "s\n" << nl;
                }
                Info<< "        " << io.name() << endl;
            }
            ++nFields;

            reconstructPointField<Type>(io)().write();
            ++nReconstructed_;
        }
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::pointFieldReconstructor::reconstructPointFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    return reconstructPointFields<Type>
    (
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    );
}


// ************************************************************************* //

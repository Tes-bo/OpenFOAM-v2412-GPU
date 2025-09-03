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

#include "fvFieldReconstructor.H"
#include "Time.H"
#include "PtrList.H"
#include "emptyFvPatch.H"
#include "fvPatchFields.H"
#include "fvsPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructField
(
    const IOobject& fieldObject,
    const PtrList<DimensionedField<Type, volMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    forAll(procMeshes_, proci)
    {
        const DimensionedField<Type, volMesh>& procField = procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.field(),
            cellProcAddressing_[proci]
        );
    }

    auto tfield = tmp<DimensionedField<Type, volMesh>>::New
    (
        fieldObject,
        mesh_,
        procFields[0].dimensions(),
        internalField
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructField
(
    const IOobject& fieldObject,
    const PtrList<GeometricField<Type, fvPatchField, volMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    // Create the patch fields
    PtrList<fvPatchField<Type>> patchFields(mesh_.boundary().size());

    forAll(procFields, proci)
    {
        const auto& procField = procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            cellProcAddressing_[proci]
        );

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[proci], patchi)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // Get addressing slice for this patch
            const labelList::subList cp =
                procField.mesh().boundary()[patchi].patchSlice
                (
                    faceProcAddressing_[proci]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields.set(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        fvPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            fvPatchField<Type>::Internal::null(),
                            fvPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                const label curPatchStart =
                    mesh_.boundaryMesh()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, facei)
                {
                    // Check
                    if (cp[facei] <= 0)
                    {
                        FatalErrorInFunction
                            << "Processor " << proci
                            << " patch "
                            << procField.mesh().boundary()[patchi].name()
                            << " face " << facei
                            << " originates from reversed face since "
                            << cp[facei]
                            << exit(FatalError);
                    }

                    // Subtract one to take into account offsets for
                    // face direction.
                    reverseAddressing[facei] = cp[facei] - 1 - curPatchStart;
                }


                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchi];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, facei)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
                    label curF = cp[facei] - 1;

                    // Is the face on the boundary?
                    if (curF >= mesh_.nInternalFaces())
                    {
                        label curBPatch = mesh_.boundaryMesh().whichPatch(curF);

                        if (!patchFields.set(curBPatch))
                        {
                            patchFields.set
                            (
                                curBPatch,
                                fvPatchField<Type>::New
                                (
                                    mesh_.boundary()[curBPatch].type(),
                                    mesh_.boundary()[curBPatch],
                                    fvPatchField<Type>::Internal::null()
                                )
                            );
                        }

                        // add the face
                        label curPatchFace =
                            mesh_.boundaryMesh()
                                [curBPatch].whichFace(curF);

                        patchFields[curBPatch][curPatchFace] =
                            curProcPatch[facei];
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchi)
    {
        // add empty patches
        if
        (
            isType<emptyFvPatch>(mesh_.boundary()[patchi])
         && !patchFields.set(patchi)
        )
        {
            patchFields.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    fvPatchFieldBase::emptyType(),
                    mesh_.boundary()[patchi],
                    fvPatchField<Type>::Internal::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    auto tfield = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        fieldObject,
        mesh_,
        procFields[0].dimensions(),
        internalField,
        patchFields
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldReconstructor::reconstructField
(
    const IOobject& fieldObject,
    const PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(mesh_.nInternalFaces());

    // Create the patch fields
    PtrList<fvsPatchField<Type>> patchFields(mesh_.boundary().size());


    forAll(procMeshes_, proci)
    {
        const auto& procField = procFields[proci];

        // Set the face values in the reconstructed field

        // It is necessary to create a copy of the addressing array to
        // take care of the face direction offset trick.
        //
        {
            const labelList& faceMap = faceProcAddressing_[proci];

            // Correctly oriented copy of internal field
            Field<Type> procInternalField(procField.primitiveField());
            // Addressing into original field
            labelList curAddr(procInternalField.size());

            forAll(procInternalField, addrI)
            {
                curAddr[addrI] = mag(faceMap[addrI])-1;
                if (faceMap[addrI] < 0)
                {
                    procInternalField[addrI] = -procInternalField[addrI];
                }
            }

            // Map
            internalField.rmap(procInternalField, curAddr);
        }

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[proci], patchi)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // Get addressing slice for this patch
            const labelList::subList cp =
                procMeshes_[proci].boundary()[patchi].patchSlice
                (
                    faceProcAddressing_[proci]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields.set(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        fvsPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            fvsPatchField<Type>::Internal::null(),
                            fvPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                const label curPatchStart =
                    mesh_.boundaryMesh()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, facei)
                {
                    // Subtract one to take into account offsets for
                    // face direction.
                    reverseAddressing[facei] = cp[facei] - 1 - curPatchStart;
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchi];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, facei)
                {
                    label curF = cp[facei] - 1;

                    // Is the face turned the right side round
                    if (curF >= 0)
                    {
                        // Is the face on the boundary?
                        if (curF >= mesh_.nInternalFaces())
                        {
                            label curBPatch =
                                mesh_.boundaryMesh().whichPatch(curF);

                            if (!patchFields.set(curBPatch))
                            {
                                patchFields.set
                                (
                                    curBPatch,
                                    fvsPatchField<Type>::New
                                    (
                                        mesh_.boundary()[curBPatch].type(),
                                        mesh_.boundary()[curBPatch],
                                        fvsPatchField<Type>::Internal::null()
                                    )
                                );
                            }

                            // add the face
                            label curPatchFace =
                                mesh_.boundaryMesh()
                                [curBPatch].whichFace(curF);

                            patchFields[curBPatch][curPatchFace] =
                                curProcPatch[facei];
                        }
                        else
                        {
                            // Internal face
                            internalField[curF] = curProcPatch[facei];
                        }
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchi)
    {
        // add empty patches
        if
        (
            isType<emptyFvPatch>(mesh_.boundary()[patchi])
         && !patchFields.set(patchi)
        )
        {
            patchFields.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    fvsPatchFieldBase::emptyType(),
                    mesh_.boundary()[patchi],
                    fvsPatchField<Type>::Internal::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    auto tfield = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
    (
        fieldObject,
        mesh_,
        procFields[0].dimensions(),
        internalField,
        patchFields
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructInternalField
(
    const IOobject& fieldObject
) const
{
    // Read the field for all the processors
    PtrList<DimensionedField<Type, volMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new DimensionedField<Type, volMesh>
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
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructVolumeField
(
    const IOobject& fieldObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, fvPatchField, volMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, fvPatchField, volMesh>
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
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldReconstructor::reconstructSurfaceField
(
    const IOobject& fieldObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, fvsPatchField, surfaceMesh>
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
Foam::label Foam::fvFieldReconstructor::reconstructInternalFields
(
    const UPtrList<const IOobject>& fieldObjects
)
{
    typedef DimensionedField<Type, volMesh> fieldType;

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

            reconstructInternalField<Type>(io)().write();
            ++nReconstructed_;
        }
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::fvFieldReconstructor::reconstructVolumeFields
(
    const UPtrList<const IOobject>& fieldObjects
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

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

            reconstructVolumeField<Type>(io)().write();
            ++nReconstructed_;
        }
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::fvFieldReconstructor::reconstructSurfaceFields
(
    const UPtrList<const IOobject>& fieldObjects
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

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


            reconstructSurfaceField<Type>(io)().write();
            ++nReconstructed_;
        }
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::fvFieldReconstructor::reconstructInternalFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef DimensionedField<Type, volMesh> fieldType;

    return reconstructInternalFields<Type>
    (
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    );
}


template<class Type>
Foam::label Foam::fvFieldReconstructor::reconstructVolumeFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    return reconstructVolumeFields<Type>
    (
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    );
}


template<class Type>
Foam::label Foam::fvFieldReconstructor::reconstructSurfaceFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    return reconstructSurfaceFields<Type>
    (
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    );
}


// ************************************************************************* //

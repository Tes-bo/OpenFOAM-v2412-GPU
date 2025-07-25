/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "Time.H"
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "distributedFieldMapper.H"
#include "distributedPointPatchFieldMapper.H"
#include "emptyPointPatch.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::parPointFieldDistributor::distributeField
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld
) const
{
    if (!tgtMeshRef_ || !distMapRef_)
    {
        FatalErrorInFunction
            << "Cannot map field without target mesh and/or distribution!"
            << abort(FatalError);
    }
    if (!hasPatchPointMaps())
    {
        const_cast<parPointFieldDistributor&>(*this).createPatchPointMaps();
    }

    const auto& tgtMesh = tgtMeshRef_();
    const auto& distMap = distMapRef_();

    // Create internalField by remote mapping
    distributedFieldMapper mapper
    (
        labelUList::null(),
        distMap.pointMap()
    );

    DimensionedField<Type, pointMesh> internalField
    (
        IOobject
        (
            fld.name(),
            tgtMesh.time().timeName(),
            fld.local(),
            tgtMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tgtMesh,
        fld.dimensions(),
        Field<Type>(fld.internalField(), mapper)
    );

    internalField.oriented() = fld.oriented();


    // Create patchFields by remote mapping

    PtrList<pointPatchField<Type>> newPatchFields(tgtMesh.boundary().size());

    const auto& bfld = fld.boundaryField();

    forAll(bfld, patchi)
    {
        if (patchPointMaps_.set(patchi))
        {
            // Clone local patch field

            const distributedPointPatchFieldMapper mapper
            (
                labelUList::null(),
                patchPointMaps_[patchi]
            );

            // Map into local copy
            newPatchFields.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    bfld[patchi],
                    tgtMesh.boundary()[patchi],    // pointPatch
                    pointPatchField<Type>::Internal::null(),
                    mapper
                )
            );

            // Note: below alternative, 'clone' method will not work since
            // there is no clone to reset both internalField reference and
            // patch reference. TBD.
            //newPatchFields.set
            //(
            //    patchi,
            //    bfld[patchi].clone
            //    (
            //        tgtMesh.boundary()[patchi],
            //        pointPatchField<Type>::Internal::null(),
            //        mapper
            //    )
            //);
        }
    }

    // Add some empty patchFields on remaining patches (this also handles
    // e.g. processorPatchFields or any other constraint type patches)
    forAll(newPatchFields, patchi)
    {
        if (!newPatchFields.set(patchi))
        {
            newPatchFields.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    pointPatchFieldBase::emptyType(),
                    tgtMesh.boundary()[patchi],
                    pointPatchField<Type>::Internal::null()
                )
            );
        }
    }

    return
        tmp<GeometricField<Type, pointPatchField, pointMesh>>::New
        (
            std::move(internalField),
            newPatchFields
        );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::parPointFieldDistributor::distributePointField
(
    const IOobject& fieldObject
) const
{
    // Read field
    GeometricField<Type, pointPatchField, pointMesh> fld
    (
        fieldObject,
        srcMesh_
    );

    // Redistribute
    return distributeField(fld);
}


template<class Type>
Foam::label Foam::parPointFieldDistributor::distributePointFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    const UPtrList<const IOobject> fieldObjects
    (
        selectedFields.empty()
      ? objects.csorted<fieldType>()
      : objects.csorted<fieldType>(selectedFields)
    );

    label nFields = 0;
    for (const IOobject& io : fieldObjects)
    {
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Reconstructing "
                    << fieldType::typeName << "s\n" << nl;
            }
            Info<< "        " << io.name() << nl;
        }
        ++nFields;

        tmp<fieldType> tfld(distributePointField<Type>(io));


        if (isWriteProc_.good())
        {
            if (isWriteProc_)
            {
                tfld().write();
            }
        }
        else if (writeHandler_ && writeHandler_->good())
        {
            auto oldHandler = fileOperation::fileHandler(writeHandler_);
            const label oldComm = UPstream::commWorld(fileHandler().comm());

            tfld().write();

            writeHandler_ = fileOperation::fileHandler(oldHandler);
            UPstream::commWorld(oldComm);
        }
    }

    if (nFields && verbose_) Info<< endl;
    return nFields;
}


template<class Type>
void Foam::parPointFieldDistributor::distributeAndStore
(
    const PtrList<GeometricField<Type, pointPatchField, pointMesh>>& fields
) const
{
    for (const auto& fld : fields)
    {
        // Distribute and store
        auto tfld = distributeField(fld);

        tfld.ref().writeOpt(IOobject::AUTO_WRITE);

        tfld().mesh().thisDb().store(tfld);
    }
}


// ************************************************************************* //

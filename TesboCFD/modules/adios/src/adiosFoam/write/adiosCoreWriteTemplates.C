/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "adiosCoreWrite.H"
#include "ListStream.H"
#include "areaFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFaPatch.H"
#include "processorFvPatch.H"
#include "processorPointPatch.H"

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

template<class FvBoundaryMesh>
Foam::label Foam::adiosFoam::adiosCoreWrite::nNonProcessor
(
    const FvBoundaryMesh& bmesh
)
{
    label nonProc = 0;

    for (const auto& p : bmesh)
    {
        if
        (
            isA<processorFaPatch>(p)
         || isA<processorFvPatch>(p)
         || isA<processorPointPatch>(p)
        )
        {
            break;
        }

        ++nonProc;
    }

    return nonProc;
}


template<class Type>
adios2::Variable<typename Foam::pTraits<Type>::cmptType>
Foam::adiosFoam::adiosCoreWrite::definePutVariable
(
    const fileName& name,
    const label local,
    const label global,
    const label offset
)
{
    using cmptType = typename pTraits<Type>::cmptType;

    adios2::Dims globals
    {
        static_cast<adios2::Dims::value_type>(global),
        static_cast<adios2::Dims::value_type>(pTraits<Type>::nComponents)
    };

    adios2::Dims offsets
    {
        static_cast<adios2::Dims::value_type>(offset)
    };

    adios2::Dims dims
    {
        static_cast<adios2::Dims::value_type>(local),
        static_cast<adios2::Dims::value_type>(pTraits<Type>::nComponents)
    };

    if (pTraits<Type>::nComponents <= 1)
    {
        dims.resize(1);
        globals.resize(1);
    }

    return
        writeIOPtr_->DefineVariable<cmptType>
        (
            name,       // shape = global 1D array of this info
            globals,    // start (offset) = for this rank
            offsets,    // single value per rank
            dims
        );
}


template<class Type>
bool Foam::adiosFoam::adiosCoreWrite::putListVariable
(
    const fileName& name,
    const UList<Type>& list
)
{
    using cmptType = typename pTraits<Type>::cmptType;

    // Type is primitive (== cmptType) or VectorSpace (contiguous)
    // containers of cmptType.
    // So it is permissible (but ugly), to use reinterpret_cast.
    // However, if the list is empty its underlying data pointer is nullptr!!

    const cmptType* data =
    (
        list.empty()
      ? &(pTraits<cmptType>::zero)
      : reinterpret_cast<const cmptType*>(list.cdata())
    );

    adios2::Dims dims
    {
        static_cast<adios2::Dims::value_type>(list.size()),
        static_cast<adios2::Dims::value_type>(pTraits<Type>::nComponents)
    };

    // Skip this and just leave nComponents == 1 in the definition too?
    if (pTraits<Type>::nComponents <= 1)
    {
        dims.resize(1);
    }

    auto var =
        writeIOPtr_->DefineVariable<cmptType>
        (
            name,
            {},     // global
            {},     // offsets
            dims    // local dims
        );

    if (writeFilePtr_)
    {
#if 0
        // DefineOperator("zfp", ...);
        if (std::is_floating_point<cmptType>())
        {
            // floating point data:
            // adios_set_transform(varid, "zfp:accuracy=0.0001");
        }
        else
        {
            // int data:  adios_set_transform(varid, "zlib");
            // adios_set_transform(varid, "bzip2");
        }
#endif

        writeFilePtr_->Put<cmptType>(var, data, adios2::Mode::Sync);
        return true;
    }

    WarningInFunction
        << "Attempting to write adios variable "
        << name << " without an open adios file"
        << endl;

    return false;
}


template<class Type>
bool Foam::adiosFoam::adiosCoreWrite::putStreamVariable
(
    const fileName& name,
    const Type& obj
)
{
    // Use transfer buffer
    auto& buf = transferBuffer();

    OListStream os(std::move(buf), IOstream::BINARY);
    os.setBlockSize(4096);
    os.indentSize() = 0;  // avoid some spaces
    os  << obj;
    os.swap(buf);

    // Write as byte-stream
    auto var =
        writeIOPtr_->DefineVariable<char>
        (
            name,
            {},
            {},
            { static_cast<adios2::Dims::value_type>(buf.size()) }
        );

    writeFilePtr_->Put<char>(var, buf.cdata(), adios2::Mode::Sync);

    return true;
}


template<class Type>
Foam::wordList Foam::adiosFoam::adiosCoreWrite::writeCloudFields
(
    const fileName& cloudVarName,
    const objectRegistry& obrTmp
)
{
    //
    // Cloud fields
    //

    const word& clsName = IOField<Type>::typeName;

    // Fields are not always on all processors (eg, multi-component parcels).
    // Thus need to resolve names between all processors.

    wordList fieldNames(obrTmp.names<IOField<Type>>());
    #if (OPENFOAM < 2212)
    Pstream::combineGather(fieldNames, ListOps::uniqueEqOp<word>());
    Pstream::combineScatter(fieldNames);
    #else
    Pstream::combineReduce(fieldNames, ListOps::uniqueEqOp<word>());
    #endif
    Foam::sort(fieldNames);  // Consistent order

    for (const word& fieldName : fieldNames)
    {
        const List<Type>* fldPtr = obrTmp.findObject<IOField<Type>>(fieldName);
        const List<Type>& values = (fldPtr ? *fldPtr : List<Type>::null());

        putListVariable(cloudVarName/fieldName, values);

        if (UPstream::master())
        {
            putAttribute(cloudVarName/fieldName/"class", clsName);
        }
    }

    return fieldNames;
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::adiosFoam::adiosCoreWrite::writeGeometricField
(
    const GeometricField<Type, PatchField, GeoMesh>* pfield
)
{
    if (pfield)
    {
        return writeGeometricField(*pfield);
    }

    return false;
}


// dimensions, internalField, boundaryField
template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::adiosFoam::adiosCoreWrite::writeGeometricField
(
    const GeometricField<Type, PatchField, GeoMesh>& field
)
{
    const fileName varName = adiosFoam::fieldPath
    (
        field.db().name(), // = region name
        field.name()       // = field name
    );

    const auto& bfield = field.boundaryField();
    const auto& bmesh  = field.mesh().boundary();
    const label nPatches = nNonProcessor(bmesh);

    if (UPstream::master())
    {
        // Attributes are global - only passed via the master
        putAttribute(varName / "class", field.type());
        putDimensionsAttribute(varName, field.dimensions());
    }

    // "internalField" as primitive field
    writeField(varName, field.internalField());

    // Independent of storage, provide a quick lookup of field patch types
    if (UPstream::master())
    {
        // Attributes are global - only passed via the master
        // - Only consider 'real' (non-processor) patches

        stringList pTypes(nPatches);
        for (label patchi=0; patchi < nPatches; ++patchi)
        {
            pTypes[patchi] = bfield[patchi].type();
        }

        putListAttribute(varName / "patch-types", pTypes);
    }

    // "boundaryField" as byte stream, but need adapter for dictionary format
    adiosFoam::BoundaryWriter<Type, PatchField, GeoMesh> adapter(bfield);

    putStreamVariable(varName / "boundaryField", adapter);

    // Also write the corresponding "_0" field if it exists
    if (field.nOldTimes() > 1)
    {
        writeGeometricField(field.oldTime());
    }

    return true;
}


template<class Type, class GeoMesh>
bool Foam::adiosFoam::adiosCoreWrite::writeInternalField
(
    const DimensionedField<Type, GeoMesh>* pfield
)
{
    if (pfield)
    {
        return writeInternalField(*pfield);
    }

    return false;
}


template<class Type, class GeoMesh>
bool Foam::adiosFoam::adiosCoreWrite::writeInternalField
(
    const DimensionedField<Type, GeoMesh>& field
)
{
    const fileName varName
    (
        adiosFoam::fieldPath
        (
            field.db().name(), // = region name
            field.name()       // = field name
        )
    );

    if (UPstream::master())
    {
        // Attributes are global - only passed via the master
        putAttribute(varName / "class", field.type());
        putDimensionsAttribute(varName, field.dimensions());
    }

    return writeField(varName, field);
}


template<class Type, class GeoMesh>
bool Foam::adiosFoam::adiosCoreWrite::writeField
(
    const fileName& varName,
    const word& local,
    const DimensionedField<Type, GeoMesh>& field
)
{
    return putListVariable(varName/local, field.field());
}


template<class Type, class GeoMesh>
bool Foam::adiosFoam::adiosCoreWrite::writeField
(
    const fileName& varName,
    const DimensionedField<Type, GeoMesh>& field
)
{
    return writeField(varName, word::null, field);
}


template<class ZoneMeshType>
void Foam::adiosFoam::adiosCoreWrite::writeMeshZones
(
    const fileName& zonePath,
    const ZoneMeshType& zones
)
{
    for (const zone& zn : zones)
    {
        putListVariable(zonePath/zn.name(), zn);
    }
}


// ************************************************************************* //

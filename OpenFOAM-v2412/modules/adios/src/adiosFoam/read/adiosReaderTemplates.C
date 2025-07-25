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

#include "adiosCore.H"
#include "adiosReader.H"
#include "ListStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class SourceType>
Foam::label Foam::adiosFoam::adiosReader::readCopyList
(
    List<Type>& list,
    const fileName& varName,
    const bool allowResize
)
{
    typedef typename pTraits<Type>::cmptType cmptType;

    adios2::Variable<SourceType> var =
        readIOPtr_->InquireVariable<SourceType>(varName);

    if (!var)
    {
        const std::string varTypeName(inquireVariableType(varName));

        if (varTypeName.empty())
        {
            FatalErrorInFunction
                << "Missing adios variable: " << varName << nl
                << exit(FatalIOError);
        }
        else
        {
            FatalErrorInFunction
                << "Type mismatch for adios variable: " << varName << nl
                << "Expecting " << adios2::GetType<SourceType>()
                << " but found " << varTypeName << nl
                << exit(FatalIOError);
        }

        return -1;
    }


    var.SetBlockSelection(UPstream::myProcNo());
    const size_t nElem = var.SelectionSize();

    if (allowResize)
    {
        list.resize(nElem / pTraits<Type>::nComponents);
    }
    else
    {
        // TODO? check sizes - Check target size?
        // (list.size() * pTraits<Type>::nComponents) == nElem
    }

    if (std::is_same<cmptType, SourceType>::value)
    {
        readFilePtr_->Get<SourceType>
        (
            var,
            reinterpret_cast<SourceType*>(list.data()),
            adios2::Mode::Sync
        );
    }
    else
    {
        // Use intermediate byte-buffer
        auto& byteBuffer = adiosCore::transferBuffer();

        byteBuffer.resize(nElem*var.Sizeof());

        UList<SourceType> inputBuffer
        (
            reinterpret_cast<SourceType*>(byteBuffer.data()),
            nElem
        );

        readFilePtr_->Get<SourceType>
        (
            var,
            inputBuffer.data(),
            adios2::Mode::Sync
        );


        // Component-wise copy. Cannot use std::copy

        auto srcIter = inputBuffer.cbegin();
        for (Type& item : list)
        {
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
            {
                setComponent(item, cmpt) = *srcIter;
                ++srcIter;
            }
        }
    }

    return (nElem / pTraits<Type>::nComponents);
}


template<class Type>
Foam::label Foam::adiosFoam::adiosReader::readList
(
    List<Type>& list,
    const fileName& varName,
    const bool allowResize
)
{
    typedef typename pTraits<Type>::cmptType SourceType;

    adios2::Variable<SourceType> var =
        readIOPtr_->InquireVariable<SourceType>(varName);

    if (!var)
    {
        const std::string varTypeName(inquireVariableType(varName));

        if (varTypeName.empty())
        {
            FatalErrorInFunction
                << "Missing adios variable: " << varName << nl
                << exit(FatalIOError);
        }
        else
        {
            FatalErrorInFunction
                << "Type mismatch for adios variable: " << varName << nl
                << "Expecting " << adios2::GetType<SourceType>()
                << " but found " << varTypeName << nl
                << exit(FatalIOError);
        }

        return -1;
    }


    return readCopyList<Type, SourceType>(list, varName, allowResize);
}


template<class Type>
Foam::label Foam::adiosFoam::adiosReader::readList
(
    List<Type>& list,
    const fileName& varName,
    const std::string& varTypeName,
    const bool allowResize
)
{
    typedef typename pTraits<Type>::cmptType cmptType;

    if (varTypeName.empty())
    {
        FatalErrorInFunction
            << "Missing adios variable: " << varName << nl
            << exit(FatalIOError);

        return -1;
    }

    label sizeRead = -1;

    if (adios2::GetType<cmptType>() == varTypeName)
    {
        // Assume the best. Need to verify later!
        sizeRead = readCopyList<Type, cmptType>(list, varName, allowResize);
    }
    else if (adios2::GetType<float>() == varTypeName)
    {
        sizeRead = readCopyList<Type, float>(list, varName, allowResize);
    }
    else if (adios2::GetType<double>() == varTypeName)
    {
        sizeRead = readCopyList<Type, double>(list, varName, allowResize);
    }
    else if (adios2::GetType<int32_t>() == varTypeName)
    {
        sizeRead = readCopyList<Type, int32_t>(list, varName, allowResize);
    }
    else if (adios2::GetType<int64_t>() == varTypeName)
    {
        sizeRead = readCopyList<Type, int64_t>(list, varName, allowResize);
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported adios type: " << varTypeName << nl
            << exit(FatalIOError);
    }

    return sizeRead;
}


// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

template<class Type, class SourceType>
bool Foam::adiosFoam::adiosReader::readAttributeImpl
(
    const string& varName,
    Type& value,
    bool mandatory
)
{
    adios2::Attribute<SourceType> var =
        readIOPtr_->InquireAttribute<SourceType>(varName);

    if (!var)
    {
        if (mandatory)
        {
            const std::string varTypeName(inquireAttributeType(varName));

            if (varTypeName.empty())
            {
                FatalErrorInFunction
                    << "Missing adios attribute: " << varName << nl
                    << exit(FatalIOError);
            }
            else
            {
                FatalErrorInFunction
                    << "Type mismatch for adios attribute: " << varName << nl
                    << "Expecting " << adios2::GetType<SourceType>()
                    << " but found " << varTypeName << nl
                    << exit(FatalIOError);
            }
        }

        return false;
    }

    const auto data = var.Data();

    if (data.size() != 1)
    {
        FatalErrorInFunction
            << "Incorrect number of elements for attribute: " << varName << nl
            << "  expecting 1 " << adios2::GetType<SourceType>()
            << ", found " << data.size() << nl
            << exit(FatalIOError);
    }

    value = data.front();

    return true;
}


template<class Type, class SourceType>
bool Foam::adiosFoam::adiosReader::readVariableImpl
(
    const string& varName,
    Type& value,
    bool mandatory
)
{
    // we probably don't distinguish between signed/unsigned
    // mostly just use signed anyhow

    adios2::Variable<SourceType> var =
        readIOPtr_->InquireVariable<SourceType>(varName);

    if (!var)
    {
        const std::string varTypeName(inquireAttributeType(varName));

        if (varTypeName.empty())
        {
            if (mandatory)
            {
                FatalErrorInFunction
                    << "Missing adios variable: " << varName << nl
                    << exit(FatalIOError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Type mismatch for adios variable: " << varName << nl
                << "Expecting " << adios2::GetType<SourceType>()
                << " but found " << varTypeName << nl
                << exit(FatalIOError);
        }

        return false;
    }

    const auto nElem = var.SelectionSize();

    if (nElem != 1)
    {
        FatalErrorInFunction
            << "wrong number of elements for variable: " << varName << nl
            << "  expecting 1 " << adios2::GetType<SourceType>()
            << ", found " << nElem << nl
            << exit(FatalIOError);
    }


    // OK for single value, but need std::copy and rewriter for dissimilar
    SourceType dummy;
    readFilePtr_->Get<SourceType>(var, dummy, adios2::Mode::Sync);

    value = dummy;

    return true;
}


template<class Type, class SourceType>
bool Foam::adiosFoam::adiosReader::readListAttributeImpl
(
    const string& varName,
    List<Type>& list,
    bool mandatory
)
{
    adios2::Attribute<SourceType> attr =
        readIOPtr_->InquireAttribute<SourceType>(varName);

    if (!attr)
    {
        const std::string varTypeName(inquireAttributeType(varName));

        if (varTypeName.empty())
        {
            if (mandatory)
            {
                FatalErrorInFunction
                    << "Missing adios variable: " << varName << nl
                    << exit(FatalIOError);
            }
        }
        else
        {
            FatalErrorInFunction
                << "Type mismatch for adios variable: " << varName << nl
                << "Expecting " << adios2::GetType<SourceType>()
                << " but found " << varTypeName << nl
                << exit(FatalIOError);
        }

        return false;
    }

    const auto data = attr.Data();

    list.resize(data.size());

    std::copy(data.begin(), data.end(), list.begin());

    return true;
}


template<class StringType>
bool Foam::adiosFoam::adiosReader::readStringListAttributeIfPresent
(
    const string& name,
    List<StringType>& list
)
{
    return readListAttributeImpl<StringType, std::string>(name, list, false);
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::adiosFoam::adiosReader::readGeometricField
(
    GeometricField<Type, PatchField, GeoMesh>* pfield,
    const fileName& varName
)
{
    if (pfield)
    {
        return readGeometricField(*pfield, varName) >= 0;
    }

    return false;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::label Foam::adiosFoam::adiosReader::readGeometricField
(
    GeometricField<Type, PatchField, GeoMesh>& field,
    const fileName& varName
)
{
    label sizeRead = readField(field.ref(), varName);

    if (adiosCore::debug & 2)
    {
        Pout<<"Read " << varName << " was " << sizeRead << nl;
    }

    size_t nread = 0;
    const fileName boundaryFieldName = varName/"boundaryField";

    adios2::Variable<char> var =
        readIOPtr_->InquireVariable<char>(boundaryFieldName);

    if (var)
    {
        var.SetBlockSelection(UPstream::myProcNo());
        nread = (var.SelectionSize()*var.Sizeof());
    }


    // Use intermediate byte-buffer
    auto& byteBuffer = adiosCore::transferBuffer();


    if (nread)
    {
        byteBuffer.resize(nread);

        readFilePtr_->Get<char>(var, byteBuffer.data(), adios2::Mode::Sync);

        // Read fields via dictionary
        UIListStream is(byteBuffer, IOstream::BINARY);
        dictionary dict(is);

        if (adiosCore::debug & 2)
        {
            Pout<<"dictionary: " << field.name() << " with "
                << dict.toc() << nl;
        }

        field.boundaryFieldRef().readField(field, dict);
    }


    // Could/should also check nread (boundaryField)

    const fileName field0Name = varName + "_0";
    if (field.nOldTimes() > 1 && hasVariable<>(field0Name))
    {
        Info<<"should load " << field.oldTime().name()
            << " from " << field0Name
            << " since we have " << field.nOldTimes() << nl;

#if 0
        // This needs more work to avoid creating time references for the
        // "_0" field too. Probably have to bypass the standard mechanism
        // and do the transfer directly and update the eventNo afterwards.
        readGeometricField(field.oldTime(), field0Name);
#endif
    }

    return sizeRead;
}


template<class Type, class GeoMesh>
bool Foam::adiosFoam::adiosReader::readField
(
    DimensionedField<Type, GeoMesh>* pfield,
    const fileName& varName
)
{
    if (pfield)
    {
        return readField(*pfield, varName) >= 0;
    }

    return false;
}


template<class Type, class GeoMesh>
Foam::label Foam::adiosFoam::adiosReader::readField
(
    DimensionedField<Type, GeoMesh>& field,
    const fileName& varName
)
{
    typedef typename pTraits<Type>::cmptType cmptType;

    const std::string varTypeName(inquireVariableType(varName));

    if (varTypeName.empty())
    {
        FatalErrorInFunction
            << "Missing adios variable: " << varName << nl
            << exit(FatalIOError);

        return -1;
    }

    Field<Type>& list = field.field();

    label sizeRead = -1;

    if (adios2::GetType<cmptType>() == varTypeName)
    {
        sizeRead = readCopyList<Type, cmptType>(list, varName);
    }
    else if (adios2::GetType<float>() == varTypeName)
    {
        sizeRead = readCopyList<Type, float>(list, varName);
    }
    else if (adios2::GetType<double>() == varTypeName)
    {
        sizeRead = readCopyList<Type, double>(list, varName);
    }
    else if (adios2::GetType<int32_t>() == varTypeName)
    {
        sizeRead = readCopyList<Type, int32_t>(list, varName);
    }
    else if (adios2::GetType<int64_t>() == varTypeName)
    {
        sizeRead = readCopyList<Type, int64_t>(list, varName);
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported adios type: " << varTypeName << nl
            << exit(FatalIOError);
    }

    return sizeRead;
}


template<class Type>
bool Foam::adiosFoam::adiosReader::readIntoRegistry
(
    objectRegistry& obr,
    const fileName& varName,
    word clsName
)
{
    if (clsName.empty())
    {
        clsName = getStringAttribute(varName/"class");
    }

    if (IOField<Type>::typeName != clsName)
    {
        return false;
    }

    auto *ptr = new IOField<Type>
    (
        IOobject
        (
            varName.name(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true  // Register
        )
    );
    ptr->store();


    const std::string varTypeName(inquireVariableType(varName));

    // Must allow resizing, but strict
    readList<Type>(*ptr, varName, varTypeName, true);

    Info<< "Read " << clsName << " size: " <<  ptr->size()
        << " : " << varName << nl;

    return true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017, 2020 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "Cloud.H"
#include "Time.H"
#include "IOPosition.H"
#include "IOdictionary.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::word Foam::Cloud<ParticleType>::cloudPropertiesName("cloudProperties");


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::readCloudUniformProperties()
{
    IOobject dictObj
    (
        cloudPropertiesName,
        time().timeName(),
        "uniform"/cloud::prefix/name(),
        db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (dictObj.typeHeaderOk<IOdictionary>(true))
    {
        const IOdictionary uniformPropsDict(dictObj);

        // Fall back to positions mode if the entry is not present for
        // backwards compatibility
        geometryType_ =
            cloud::geometryTypeNames.getOrDefault
            (
                "geometry",
                uniformPropsDict,
                cloud::geometryType::POSITIONS
            );

        const word procName("processor" + Foam::name(Pstream::myProcNo()));

        const dictionary* dictptr = uniformPropsDict.findDict(procName);

        if (dictptr)
        {
            dictptr->readEntry("particleCount", ParticleType::particleCount_);
        }
    }
    else
    {
        ParticleType::particleCount_ = 0;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeCloudUniformProperties() const
{
    IOdictionary uniformPropsDict
    (
        IOobject
        (
            cloudPropertiesName,
            time().timeName(),
            "uniform"/cloud::prefix/name(),
            db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    labelList np(UPstream::nProcs(), Foam::zero{});
    np[UPstream::myProcNo()] = ParticleType::particleCount_;
    Pstream::allGatherList(np);

    uniformPropsDict.add
    (
        "geometry",
        cloud::geometryTypeNames[geometryType_]
    );

    forAll(np, i)
    {
        const word procName("processor" + Foam::name(i));
        uniformPropsDict.subDictOrAdd(procName).add("particleCount", np[i]);
    }

    uniformPropsDict.writeObject
    (
        IOstreamOption(IOstreamOption::ASCII, time().writeCompression()),
        true
    );
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::initCloud(const bool checkClass)
{
    readCloudUniformProperties();

    IOPosition<Cloud<ParticleType>> ioP(*this, geometryType_);

    const bool haveFile = ioP.headerOk();
    Istream& is = ioP.readStream(checkClass ? typeName : word::null, haveFile);
    if (haveFile)
    {
        ioP.readData(is, *this);
        ioP.close();
    }

    if (!haveFile && debug)
    {
        Pout<< "Not reading particle positions file: "
            << ioP.objectRelPath() << nl
            << "Assuming the initial cloud contains 0 particles." << endl;
    }

    // Always operate in coordinates mode after reading
    geometryType_ = cloud::geometryType::COORDINATES;

    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    (void)polyMesh_.tetBasePtIs();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const bool checkClass
)
:
    Cloud<ParticleType>(pMesh, Foam::zero{}, cloudName)
{
    initCloud(checkClass);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldIOobject
(
    const Cloud<ParticleType>& c,
    const IOField<DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorInFunction
            << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldFieldIOobject
(
    const Cloud<ParticleType>& c,
    const CompactIOField<Field<DataType>, DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorInFunction
            << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
template<class Type>
bool Foam::Cloud<ParticleType>::readStoreFile
(
    const IOobject& io,
    const IOobject& ioNew
) const
{
    if (io.isHeaderClass<IOField<Type>>())
    {
        IOField<Type> fld(io);
        auto* fldNewPtr = new IOField<Type>(ioNew, std::move(fld));
        return fldNewPtr->store();
    }

    return false;
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::readFromFiles
(
    objectRegistry& obr,
    const wordRes& selectFields,
    const wordRes& excludeFields
) const
{
    IOobjectList cloudObjects
    (
        *this,
        time().timeName(),
        IOobject::NO_REGISTER
    );

    const wordRes::filter pred(selectFields, excludeFields);

    forAllConstIters(cloudObjects, iter)
    {
        const IOobject& io = *(iter.val());
        const word& fldName = io.name();

        if (!pred(fldName))
        {
            continue;  // reject
        }

        IOobject ioNew
        (
            fldName,
            time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        const bool stored
        (
            readStoreFile<label>(io, ioNew)
         || readStoreFile<scalar>(io, ioNew)
         || readStoreFile<vector>(io, ioNew)
         || readStoreFile<sphericalTensor>(io, ioNew)
         || readStoreFile<symmTensor>(io, ioNew)
         || readStoreFile<tensor>(io, ioNew)
        );

        if (!stored)
        {
            DebugInfo
                << "Unhandled field:" << fldName
                << " type:" << io.headerClassName() << endl;
        }
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeFields() const
{
    ParticleType::writeFields(*this);
}


template<class ParticleType>
bool Foam::Cloud<ParticleType>::writeObject
(
    IOstreamOption streamOpt,
    const bool /* writeOnProc */
) const
{
    writeCloudUniformProperties();

    writeFields();
    return cloud::writeObject(streamOpt, (this->size() > 0));
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Cloud<ParticleType>& c)
{
    c.writeData(os);

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

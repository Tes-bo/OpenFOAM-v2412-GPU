/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "mappedPatchFieldBase.H"
#include "mappedPatchBase.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Type Foam::mappedPatchFieldBase<Type>::getAverage
(
    const dictionary& dict,
    const bool mandatory
)
{
    if (mandatory)
    {
        return dict.get<Type>("average");
    }

    return Zero;
}


template<class Type>
template<class T>
void Foam::mappedPatchFieldBase<Type>::storeField
(
    const objectRegistry& obr,
    const word& region,
    const word& patch,
    const label myComm,
    const labelListList& procToMap,
    const word& fieldName,
    const Field<T>& fld
) const
{
    // Store my data onto database

    const auto& procIDs = UPstream::procID(myComm);

    forAll(procToMap, ranki)
    {
        const labelList& map = procToMap[ranki];
        const label proci = procIDs[ranki];

        if (map.size())
        {
            const Field<T> subFld(fld, map);

            auto& subObr = const_cast<objectRegistry&>
            (
                mappedPatchBase::subRegistry
                (
                    obr,
                    mapper_.sendPath(proci)
                  / region
                  / patch
                )
            );

            if (fvPatchField<Type>::debug)
            {
                Pout<< "*** STORING :"
                    << " field:" << fieldName
                    << " values:" << flatOutput(subFld)
                    << " as:" << subObr.objectPath() << endl;
            }

            mappedPatchBase::storeField(subObr, fieldName, subFld);
        }
    }
}


template<class Type>
template<class T>
bool Foam::mappedPatchFieldBase<Type>::retrieveField
(
    const bool allowUnset,
    const objectRegistry& obr,
    const word& region,
    const word& patch,
    const label myComm,
    const labelListList& procToMap,
    const word& fieldName,
    Field<T>& fld
) const
{
    const auto& procIDs = UPstream::procID(myComm);

    bool ok = true;

    forAll(procToMap, ranki)
    {
        const labelList& map = procToMap[ranki];
        const label proci = procIDs[ranki];

        if (map.size())
        {
            auto& subObr = const_cast<objectRegistry&>
            (
                mappedPatchBase::subRegistry
                (
                    obr,
                    mapper_.receivePath(proci)
                  / region
                  / patch
                )
            );

            const IOField<T>* subFldPtr = subObr.getObjectPtr<IOField<T>>
            (
                fieldName
            );
            if (subFldPtr)
            {
                if (subFldPtr->size() != map.size())
                {
                    // This is the dummy value inserted at start-up since the
                    // map is always non-zero size (checked above)
                    //Pout<< "*** RETRIEVED DUMMY :"
                    //    << " field:" << fieldName
                    //    << " subFldPtr:" << subFldPtr->size()
                    //    << " map:" << map.size() << endl;

                    ok = false;
                }
                else
                {
                    UIndirectList<T>(fld, map) = *subFldPtr;

                    if (fvPatchField<Type>::debug)
                    {
                        Pout<< "*** RETRIEVED :"
                            << " field:" << fieldName
                            << " values:" << flatOutput(fld)
                            << " from:" << subObr.objectPath() << endl;
                    }
                }
            }
            else if (allowUnset)
            {
                if (fvPatchField<Type>::debug)
                {
                    WarningInFunction << "Not found"
                        << " field:" << fieldName
                        << " in:" << subObr.objectPath() << endl;
                }

                // Store dummy value so the database has something on it.
                // Note that size 0 should never occur naturally so we can
                // detect it if necessary.
                const Field<T> dummyFld(0);

                mappedPatchBase::storeField(subObr, fieldName, dummyFld);

                ok = false;
            }
            else
            {
                // Not found. Make it fail
                (void)subObr.lookupObject<IOField<T>>(fieldName);
                ok = false;
            }
        }
    }
    return ok;
}


template<class Type>
template<class T>
void Foam::mappedPatchFieldBase<Type>::initRetrieveField
(
    const objectRegistry& obr,
    const word& region,
    const word& patch,
    const labelListList& map,
    const word& fieldName,
    const Field<T>& fld
) const
{
    // Old code. Likely not quite correct...

    // Store my data onto database
    const label nProcs = Pstream::nProcs(0);    // comm_

    for (label domain = 0; domain < nProcs; domain++)
    {
        const labelList& constructMap = map[domain];
        if (constructMap.size())
        {
            auto& subObr = const_cast<objectRegistry&>
            (
                mappedPatchBase::subRegistry
                (
                    obr,
                    mapper_.receivePath(domain)
                  / region
                  / patch
                )
            );

            const Field<T> receiveFld(fld, constructMap);

            if (fvPatchField<Type>::debug)
            {
                Pout<< "*** STORING INITIAL :"
                    << " field:" << fieldName << " values:"
                    << flatOutput(receiveFld)
                    << " from:" << flatOutput(fld)
                    << " constructMap:" << flatOutput(constructMap)
                    << " as:" << subObr.objectPath() << endl;
            }

            mappedPatchBase::storeField(subObr, fieldName, receiveFld);
        }
    }
}


template<class Type>
template<class T>
bool Foam::mappedPatchFieldBase<Type>::storeAndRetrieveField
(
    const word& fieldName,
    const label myComm,
    const labelListList& subMap,
    const label constructSize,
    const labelListList& constructMap,
    const labelListList& address,
    const scalarListList& weights,
    Field<T>& fld
) const
{
    storeField
    (
        patchField_.internalField().time(),
        patchField_.patch().boundaryMesh().mesh().name(),
        patchField_.patch().name(),
        myComm,
        subMap,
        fieldName,
        fld
    );

    Field<T> work(constructSize);
    const bool ok = retrieveField
    (
        true,                           // allow unset
        patchField_.internalField().time(),
        mapper_.sampleRegion(),
        mapper_.samplePatch(),
        myComm,
        constructMap,
        fieldName,
        work
    );

    if (ok)
    {
        // Do interpolation

        fld.setSize(address.size());
        fld = Zero;

        const plusEqOp<T> cop;
        const multiplyWeightedOp<T, plusEqOp<T>> mop(cop);

        forAll(address, facei)
        {
            const labelList& slots = address[facei];
            const scalarList& w = weights[facei];

            forAll(slots, i)
            {
                mop(fld[facei], facei, work[slots[i]], w[i]);
            }
        }
    }
    else
    {
        // Leave fld intact
    }

    return ok;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const word& fieldName,
    const bool setAverage,
    const Type average,
    const word& interpolationScheme
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(fieldName),
    setAverage_(setAverage),
    average_(average),
    interpolationScheme_(interpolationScheme)
{}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const dictionary& dict
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_
    (
        dict.template getOrDefault<word>
        (
            "field",
            patchField_.internalField().name()
        )
    ),
    setAverage_(dict.getOrDefault("setAverage", false)),
    average_(getAverage(dict, setAverage_)),
    interpolationScheme_(interpolationCell<Type>::typeName)
{
    if
    (
        mapper_.sampleDatabase()
     && (
            mapper_.mode() != mappedPatchBase::NEARESTPATCHFACE
         && mapper_.mode() != mappedPatchBase::NEARESTPATCHFACEAMI
        )
    )
    {
        FatalErrorInFunction
            << "Mapping using the database only supported for "
            << "sampleModes "
            <<  mappedPatchBase::sampleModeNames_
                [
                    mappedPatchBase::NEARESTPATCHFACE
                ]
            << " and "
            <<  mappedPatchBase::sampleModeNames_
                [
                    mappedPatchBase::NEARESTPATCHFACEAMI
                ]
            << exit(FatalError);
    }

    if (mapper_.mode() == mappedPatchBase::NEARESTCELL)
    {
        dict.readEntry("interpolationScheme", interpolationScheme_);
    }

    // Note: in database mode derived boundary conditions need to initialise
    //       fields
}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const dictionary& dict,
    const Field<Type>& fld
)
:
    mappedPatchFieldBase<Type>::mappedPatchFieldBase(mapper, patchField, dict)
{
    if (mapper_.sampleDatabase())
    {
        if (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACE)
        {
            // Store my data on receive buffers so we have some initial data
            initRetrieveField
            (
                patchField_.internalField().time(),
                //patchField_.patch().boundaryMesh().mesh().name(),
                mapper_.sampleRegion(),
                //patchField_.patch().name(),
                mapper_.samplePatch(),
                mapper_.map().constructMap(),
                patchField_.internalField().name(),
                patchField_
            );
        }
        else if (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACEAMI)
        {
            // Depend on fall-back (sorting dummy field) in retrieveField
            // since it would be too hard to determine the field that gives
            // the wanted result after interpolation
        }
    }
}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(patchField_.internalField().name()),
    setAverage_(false),
    average_(Zero),
    interpolationScheme_(interpolationCell<Type>::typeName)
{}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchFieldBase<Type>& mapper
)
:
    mapper_(mapper.mapper_),
    patchField_(mapper.patchField_),
    fieldName_(mapper.fieldName_),
    setAverage_(mapper.setAverage_),
    average_(mapper.average_),
    interpolationScheme_(mapper.interpolationScheme_)
{}


template<class Type>
Foam::mappedPatchFieldBase<Type>::mappedPatchFieldBase
(
    const mappedPatchBase& mapper,
    const fvPatchField<Type>& patchField,
    const mappedPatchFieldBase<Type>& base
)
:
    mapper_(mapper),
    patchField_(patchField),
    fieldName_(base.fieldName_),
    setAverage_(base.setAverage_),
    average_(base.average_),
    interpolationScheme_(base.interpolationScheme_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<class Type2>
const Foam::GeometricField<Type2, Foam::fvPatchField, Foam::volMesh>&
Foam::mappedPatchFieldBase<Type>::sampleField(const word& fieldName) const
{
    typedef GeometricField<Type2, fvPatchField, volMesh> fieldType;

    if (mapper_.sameRegion())
    {
        if (fieldName == patchField_.internalField().name())
        {
            // Optimisation: bypass field lookup
            return
                dynamic_cast<const fieldType&>
                (
                    patchField_.internalField()
                );
        }
        else
        {
            const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();
            return thisMesh.template lookupObject<fieldType>(fieldName);
        }
    }

    const fvMesh& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());
    return nbrMesh.template lookupObject<fieldType>(fieldName);
}


template<class Type>
const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&
Foam::mappedPatchFieldBase<Type>::sampleField() const
{
    return sampleField<Type>(fieldName_);
}


template<class Type>
template<class T>
void Foam::mappedPatchFieldBase<Type>::distribute
(
    const word& fieldName,
    Field<T>& fld
) const
{
    if (mapper_.sampleDatabase())
    {
        const label myComm = mapper_.getCommunicator();  // Get or create

        if (mapper_.mode() != mappedPatchBase::NEARESTPATCHFACEAMI)
        {
            // Store my data on send buffers
            storeField
            (
                patchField_.internalField().time(),
                patchField_.patch().boundaryMesh().mesh().name(),
                patchField_.patch().name(),
                myComm,
                mapper_.map().subMap(),
                fieldName,
                fld
            );
            // Construct my data from receive buffers
            fld.setSize(mapper_.map().constructSize());
            retrieveField
            (
                true,                           // allow unset
                patchField_.internalField().time(),
                mapper_.sampleRegion(),
                mapper_.samplePatch(),
                myComm,
                mapper_.map().constructMap(),
                fieldName,
                fld
            );
        }
        else
        {
            const AMIPatchToPatchInterpolation& AMI = mapper_.AMI();

            // The AMI does an interpolateToSource/ToTarget. This is a
            // mapDistribute (so using subMap/constructMap) and then a
            // weighted sum. We'll store the sent data as before and
            // do the weighted summation after the retrieveField

            if (mapper_.masterWorld())
            {
                // See AMIInterpolation::interpolateToSource. Use tgtMap,
                // srcAddress, srcWeights
                storeAndRetrieveField
                (
                    fieldName,
                    myComm,
                    AMI.srcMap().subMap(),
                    AMI.tgtMap().constructSize(),
                    AMI.tgtMap().constructMap(),
                    AMI.srcAddress(),
                    AMI.srcWeights(),
                    fld
                );
            }
            else
            {
                // See AMIInterpolation::interpolateToTarget.
                // Use srcMap, tgtAddress, tgtWeights
                storeAndRetrieveField
                (
                    fieldName,
                    myComm,
                    AMI.tgtMap().subMap(),
                    AMI.srcMap().constructSize(),
                    AMI.srcMap().constructMap(),
                    AMI.tgtAddress(),
                    AMI.tgtWeights(),
                    fld
                );
            }
        }
    }
    else
    {
        mapper_.distribute(fld);
    }
}


template<class Type>
//template<class T>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchFieldBase<Type>::mappedField
(
//    const GeometricField<T, fvPatchField, volMesh>& fld
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    const fvMesh& thisMesh = patchField_.patch().boundaryMesh().mesh();

    // Result of obtaining remote values
    auto tnewValues = tmp<Field<Type>>::New();
    auto& newValues = tnewValues.ref();

    switch (mapper_.mode())
    {
        case mappedPatchBase::NEARESTCELL:
        {
            const fieldType& fld = sampleField();
            const mapDistribute& distMap = mapper_.map();

            if (interpolationScheme_ != interpolationCell<Type>::typeName)
            {
                if (!mapper_.sameWorld() || mapper_.sampleDatabase())
                {
                    FatalErrorInFunction
                        << "Interpolating cell values from different world"
                        << " or database currently not supported"
                        << exit(FatalError);
                }

                const fvMesh& nbrMesh =
                    refCast<const fvMesh>(mapper_.sampleMesh());

                // Send back sample points to the processor that holds the cell
                vectorField samples(mapper_.samplePoints());

                distMap.reverseDistribute
                (
                    (
                        mapper_.sameRegion()
                      ? thisMesh.nCells()
                      : nbrMesh.nCells()
                    ),
                    point::max,
                    samples
                );

                auto interpolator =
                    interpolation<Type>::New
                    (
                        interpolationScheme_,
                        fld
                    );

                const auto& interp = *interpolator;

                newValues.setSize(samples.size(), pTraits<Type>::max);
                forAll(samples, celli)
                {
                    if (samples[celli] != point::max)
                    {
                        newValues[celli] = interp.interpolate
                        (
                            samples[celli],
                            celli
                        );
                    }
                }
            }
            else
            {
                newValues = fld;
            }

            distribute(fieldName_, newValues);

            break;
        }
        case mappedPatchBase::NEARESTPATCHFACE:
        case mappedPatchBase::NEARESTPATCHFACEAMI:
        {
            if (mapper_.sameWorld())
            {
                const fvMesh& nbrMesh =
                    refCast<const fvMesh>(mapper_.sampleMesh());
                const fieldType& fld = sampleField();

                const label nbrPatchID =
                    nbrMesh.boundaryMesh().findPatchID(mapper_.samplePatch());

                if (nbrPatchID < 0)
                {
                    FatalErrorInFunction
                     << "Unable to find sample patch " << mapper_.samplePatch()
                     << " in region " << mapper_.sampleRegion()
                     << " for patch " << patchField_.patch().name() << nl
                     << abort(FatalError);
                }

                const auto& nbrField = fld;

                newValues = nbrField.boundaryField()[nbrPatchID];
            }
            else
            {
                // Start off from my patch values, let distribute function below
                // do all the work
                newValues = patchField_;
            }
            distribute(fieldName_, newValues);

            break;
        }
        case mappedPatchBase::NEARESTFACE:
        {
            Field<Type> allValues;
            if (mapper_.sameWorld())
            {
                const fvMesh& nbrMesh =
                    refCast<const fvMesh>(mapper_.sampleMesh());
                const fieldType& fld = sampleField();

                allValues.setSize(nbrMesh.nFaces(), Zero);

                const auto& nbrField = fld;

                for (const fvPatchField<Type>& pf : nbrField.boundaryField())
                {
                    label faceStart = pf.patch().start();

                    forAll(pf, facei)
                    {
                        allValues[faceStart++] = pf[facei];
                    }
                }
            }
            else
            {
                // Start off from my patch values. Let distribute function below
                // do all the work
                allValues.setSize(thisMesh.nFaces(), Zero);

                const fieldType& thisFld = dynamic_cast<const fieldType&>
                (
                    patchField_.internalField()
                );

                for (const fvPatchField<Type>& pf : thisFld.boundaryField())
                {
                    label faceStart = pf.patch().start();

                    forAll(pf, facei)
                    {
                        allValues[faceStart++] = pf[facei];
                    }
                }
            }

            distribute(fieldName_, allValues);
            newValues.transfer(allValues);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown sampling mode: " << mapper_.mode() << nl
                << abort(FatalError);
        }
    }

    if (setAverage_)
    {
        Type averagePsi =
            gSum(patchField_.patch().magSf()*newValues)
           /gSum(patchField_.patch().magSf());

        if (mag(averagePsi) > 0.5*mag(average_))
        {
            newValues *= mag(average_)/mag(averagePsi);
        }
        else
        {
            newValues += (average_ - averagePsi);
        }
    }

    UPstream::msgType(oldTag);  // Restore tag

    return tnewValues;
}


//template<class Type>
//Foam::tmp<Foam::Field<Type>>
//Foam::mappedPatchFieldBase<Type>::mappedField() const
//{
//    const GeometricField<Type, fvPatchField, volMesh>& fld = sampleField();
//    return mappedField<Type>(fld);
//}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedPatchFieldBase<Type>::mappedInternalField() const
{
    // Swap to obtain full local values of neighbour internal field
    auto tnbrIntFld = tmp<Field<Type>>::New();
    auto& nbrIntFld = tnbrIntFld.ref();

    if (mapper_.sameWorld())
    {
        // Same world so lookup
        const label nbrPatchID = mapper_.samplePolyPatch().index();
        const auto& nbrField = this->sampleField();
        nbrIntFld = nbrField.boundaryField()[nbrPatchID].patchInternalField();
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering
        nbrIntFld = patchField_.patchInternalField();
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    distribute(fieldName_, nbrIntFld);

    UPstream::msgType(oldTag);  // Restore tag

    return tnbrIntFld;
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::mappedPatchFieldBase<Type>::mappedWeightField() const
{
    // Swap to obtain full local values of neighbour internal field
    auto tnbrKDelta = tmp<scalarField>::New();
    auto& nbrKDelta = tnbrKDelta.ref();

    if (mapper_.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(this->mapper_.sampleMesh());
        const label nbrPatchID = mapper_.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];
        nbrKDelta = nbrPatch.deltaCoeffs();
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering
        nbrKDelta = patchField_.patch().deltaCoeffs();
    }


    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    distribute(fieldName_ + "_deltaCoeffs", nbrKDelta);

    UPstream::msgType(oldTag);  // Restore tag

    return tnbrKDelta;
}


template<class Type>
void Foam::mappedPatchFieldBase<Type>::mappedWeightField
(
    const word& fieldName,
    tmp<scalarField>& thisWeights,
    tmp<scalarField>& nbrWeights
) const
{
    thisWeights = new scalarField(patchField_.patch().deltaCoeffs());
    if (!fieldName.empty())
    {
        thisWeights.ref() *=
            patchField_.patch().template lookupPatchField<volScalarField>
            (
                fieldName
            ).patchInternalField();
    }


    // Swap to obtain full local values of neighbour internal field

    if (mapper_.sameWorld())
    {
        // Same world so lookup
        const auto& nbrMesh = refCast<const fvMesh>(mapper_.sampleMesh());
        const label nbrPatchID = mapper_.samplePolyPatch().index();
        const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

        nbrWeights = new scalarField(nbrPatch.deltaCoeffs());

        if (!fieldName.empty())
        {
            // Weightfield is volScalarField
            const auto& nbrWeightField =
                nbrMesh.template lookupObject<volScalarField>(fieldName);
            nbrWeights.ref() *=
                nbrWeightField.boundaryField()[nbrPatchID].patchInternalField();
        }
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering
        nbrWeights = new scalarField(thisWeights());
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    distribute(fieldName_ + "_weights", nbrWeights.ref());

    UPstream::msgType(oldTag);  // Restore tag
}


template<class Type>
const Foam::mappedPatchBase& Foam::mappedPatchFieldBase<Type>::mapper
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "Incorrect patch type " << p.patch().type()
            << " for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath() << nl
            << "Type should be a mappedPatch"
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}


template<class Type>
template<class T>
void Foam::mappedPatchFieldBase<Type>::initRetrieveField
(
    const word& fieldName,
    const Field<T>& fld
) const
{
    if (mapper_.sampleDatabase())
    {
        // Store my data on receive buffers (reverse of storeField;
        // i.e. retrieveField will obtain patchField)
        if (mapper_.mode() == mappedPatchBase::NEARESTPATCHFACE)
        {
            initRetrieveField
            (
                patchField_.internalField().time(),
                mapper_.sampleRegion(),
                mapper_.samplePatch(),
                mapper_.map().constructMap(),
                fieldName,
                fld
            );
        }
    }
}


template<class Type>
void Foam::mappedPatchFieldBase<Type>::write(Ostream& os) const
{
    //os.writeEntry("field", fieldName_);
    os.writeEntryIfDifferent<word>
    (
        "field",
        patchField_.internalField().name(),
        fieldName_
    );

    if (setAverage_)
    {
        os.writeEntry("setAverage", "true");
        os.writeEntry("average", average_);
    }

    if (mapper_.mode() == mappedPatchBase::NEARESTCELL)
    {
        os.writeEntry("interpolationScheme", interpolationScheme_);
    }
}


// ************************************************************************* //

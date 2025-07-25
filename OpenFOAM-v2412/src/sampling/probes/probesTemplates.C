/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "probes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOmanip.H"
#include "interpolation.H"
#include "SpanStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T>
struct isNotEqOp
{
    void operator()(T& x, const T& y) const
    {
        const T unsetVal(-VGREAT*pTraits<T>::one);

        if (x != unsetVal)
        {
            // Keep x.

            // Note: should check for y != unsetVal but multiple sample cells
            // already handled in read().
        }
        else
        {
            // x is not set. y might be.
            x = y;
        }
    }
};

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<GeoField>
Foam::probes::getOrLoadField(const word& fieldName) const
{
    tmp<GeoField> tfield;

    if (loadFromFiles_)
    {
        tfield.emplace
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobjectOption::MUST_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::NO_REGISTER
            ),
            mesh_
        );
    }
    else
    {
        tfield.cref(mesh_.cfindObject<GeoField>(fieldName));
    }

    return tfield;
}


template<class Type>
void Foam::probes::storeResults
(
    const word& fieldName,
    const Field<Type>& values
)
{
    const MinMax<Type> limits(values);
    const Type avgVal = average(values);

    this->setResult("average(" + fieldName + ")", avgVal);
    this->setResult("min(" + fieldName + ")", limits.min());
    this->setResult("max(" + fieldName + ")", limits.max());
    this->setResult("size(" + fieldName + ")", values.size());

    if (verbose_)
    {
        Info<< name() << " : " << fieldName << nl
            << "    avg: " << avgVal << nl
            << "    min: " << limits.min() << nl
            << "    max: " << limits.max() << nl << nl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::probes::writeValues
(
    const word& fieldName,
    const Field<Type>& values,
    const scalar timeValue
)
{
    if (Pstream::master())
    {
        const unsigned int width(IOstream::defaultPrecision() + 7);
        OFstream& os = *probeFilePtrs_[fieldName];

        os  << setw(width) << timeValue;

        OCharStream buf;

        forAll(values, probei)
        {
            if (includeOutOfBounds_ || processor_[probei] != -1)
            {
                buf.rewind();
                buf << values[probei];
                os  << ' ' << setw(width) << buf.str().data();
            }
        }
        os  << endl;
    }
}


template<class GeoField>
void Foam::probes::performAction
(
    const fieldGroup<GeoField>& fieldNames,
    unsigned request
)
{
    for (const word& fieldName : fieldNames)
    {
        tmp<GeoField> tfield = getOrLoadField<GeoField>(fieldName);

        if (tfield)
        {
            const auto& field = tfield();
            const scalar timeValue = field.time().timeOutputValue();

            Field<typename GeoField::value_type> values(sample(field));

            this->storeResults(fieldName, values);
            if (request & ACTION_WRITE)
            {
                this->writeValues(fieldName, values, timeValue);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample(const VolumeField<Type>& vField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    if (fixedLocations_)
    {
        autoPtr<interpolation<Type>> interpPtr
        (
            interpolation<Type>::New(samplePointScheme_, vField)
        );

        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                const vector& position = operator[](probei);

                values[probei] = interpPtr().interpolate
                (
                    position,
                    elementList_[probei],
                    -1
                );
            }
        }
    }
    else
    {
        forAll(*this, probei)
        {
            if (elementList_[probei] >= 0)
            {
                values[probei] = vField[elementList_[probei]];
            }
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample(const SurfaceField<Type>& sField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    forAll(*this, probei)
    {
        if (faceList_[probei] >= 0)
        {
            values[probei] = sField[faceList_[probei]];
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sample(const word& fieldName) const
{
    return sample(mesh_.lookupObject<VolumeField<Type>>(fieldName));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::probes::sampleSurfaceField(const word& fieldName) const
{
    return sample(mesh_.lookupObject<SurfaceField<Type>>(fieldName));
}


// ************************************************************************* //

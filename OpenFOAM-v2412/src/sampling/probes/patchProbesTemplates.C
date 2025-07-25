/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "patchProbes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::patchProbes::writeValues
(
    const word& fieldName,
    const Field<Type>& values,
    const scalar timeValue
)
{
    if (Pstream::master())
    {
        const unsigned int w = IOstream::defaultPrecision() + 7;
        OFstream& os = *probeFilePtrs_[fieldName];

        os  << setw(w) << timeValue;

        for (const auto& v : values)
        {
            os  << ' ' << setw(w) << v;
        }
        os  << endl;
    }
}


template<class GeoField>
void Foam::patchProbes::performAction
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
Foam::patchProbes::sample(const VolumeField<Type>& vField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const auto& bField = vField.boundaryField();

    forAll(*this, probei)
    {
        label facei = faceList_[probei];

        if (facei >= 0)
        {
            label patchi = patches.whichPatch(facei);
            label localFacei = patches[patchi].whichFace(facei);
            values[probei] = bField[patchi][localFacei];
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample(const SurfaceField<Type>& sField) const
{
    const Type unsetVal(-VGREAT*pTraits<Type>::one);

    auto tvalues = tmp<Field<Type>>::New(Field<Type>(this->size(), unsetVal));
    auto& values = tvalues.ref();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const auto& bField = sField.boundaryField();

    forAll(*this, probei)
    {
        label facei = faceList_[probei];

        if (facei >= 0)
        {
            label patchi = patches.whichPatch(facei);
            label localFacei = patches[patchi].whichFace(facei);
            values[probei] = bField[patchi][localFacei];
        }
    }

    Pstream::listCombineReduce(values, isNotEqOp<Type>());

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sample(const word& fieldName) const
{
    return sample(mesh_.lookupObject<VolumeField<Type>>(fieldName));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchProbes::sampleSurfaceField(const word& fieldName) const
{
    return sample(mesh_.lookupObject<SurfaceField<Type>>(fieldName));
}


// ************************************************************************* //

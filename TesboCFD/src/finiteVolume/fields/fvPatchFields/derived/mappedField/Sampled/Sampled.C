/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "fvMesh.H"
#include "volFields.H"
#include "interpolationCell.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Type Foam::PatchFunction1Types::Sampled<Type>::getAverage
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::Sampled<Type>::Sampled
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    mappedPatchBase(pp, dict),
    fieldName_(dict.get<word>("field")),
    setAverage_(dict.getOrDefault("setAverage", false)),
    average_(getAverage(dict, setAverage_)),
    interpolationScheme_(interpolationCell<Type>::typeName)
{
    if (this->mode() == mappedPatchBase::NEARESTCELL)
    {
        dict.readEntry("interpolationScheme", interpolationScheme_);
    }
}


template<class Type>
Foam::PatchFunction1Types::Sampled<Type>::Sampled
(
    const Sampled<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    mappedPatchBase(pp, rhs),
    fieldName_(rhs.fieldName_),
    setAverage_(rhs.setAverage_),
    average_(rhs.average_),
    interpolationScheme_(rhs.interpolationScheme_)
{}


template<class Type>
Foam::PatchFunction1Types::Sampled<Type>::Sampled
(
    const Sampled<Type>& rhs
)
:
    Sampled<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&
Foam::PatchFunction1Types::Sampled<Type>::sampleField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (this->sameRegion())
    {
        const polyMesh& thisMesh =
            this->mappedPatchBase::patch_.boundaryMesh().mesh();
        return thisMesh.template lookupObject<fieldType>(fieldName_);
    }
    else
    {
        const fvMesh& nbrMesh = refCast<const fvMesh>(this->sampleMesh());
        return nbrMesh.template lookupObject<fieldType>(fieldName_);
    }
}


template<class Type>
bool Foam::PatchFunction1Types::Sampled<Type>::haveSampleField() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (this->sameRegion())
    {
        const polyMesh& thisMesh =
            this->mappedPatchBase::patch_.boundaryMesh().mesh();
        return thisMesh.template foundObject<fieldType>(fieldName_);
    }
    else
    {
        const fvMesh& nbrMesh = refCast<const fvMesh>(this->sampleMesh());
        return nbrMesh.template foundObject<fieldType>(fieldName_);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::Sampled<Type>::value
(
    const scalar x
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    const fvMesh& thisMesh = refCast<const fvMesh>
    (
        this->mappedPatchBase::patch_.boundaryMesh().mesh()
    );
    const fvMesh& nbrMesh = refCast<const fvMesh>(this->sampleMesh());


    // Result of obtaining remote values
    auto tnewValues = tmp<Field<Type>>::New();
    auto& newValues = tnewValues.ref();

    if (!haveSampleField())
    {
        UPstream::msgType(oldTag);  // Restore tag
        newValues.resize_nocopy(this->mappedPatchBase::patch_.size());
        newValues = Zero;
        return this->transform(tnewValues);
    }

    switch (this->mode())
    {
        case mappedPatchBase::NEARESTCELL:
        {
            const mapDistribute& distMap = this->map();

            if (interpolationScheme_ != interpolationCell<Type>::typeName)
            {
                // Send back sample points to the processor that holds the cell
                vectorField samples(this->samplePoints());
                distMap.reverseDistribute
                (
                    (
                        this->sameRegion()
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
                        sampleField()
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
                newValues = sampleField();
            }
            distMap.distribute(newValues);

            break;
        }
        case mappedPatchBase::NEARESTPATCHFACE:
        case mappedPatchBase::NEARESTPATCHFACEAMI:
        {
            const label nbrPatchID =
                nbrMesh.boundaryMesh().findPatchID(this->samplePatch());

            if (nbrPatchID < 0)
            {
                FatalErrorInFunction
                 << "Unable to find sample patch " << this->samplePatch()
                 << " in region " << this->sampleRegion()
                 << " for patch " << this->mappedPatchBase::patch_.name() << nl
                 << abort(FatalError);
            }

            const fieldType& nbrField = sampleField();

            newValues = nbrField.boundaryField()[nbrPatchID];
            this->distribute(newValues);

            break;
        }
        case mappedPatchBase::NEARESTFACE:
        {
            Field<Type> allValues(nbrMesh.nFaces(), Zero);

            const fieldType& nbrField = sampleField();

            for (const fvPatchField<Type>& pf : nbrField.boundaryField())
            {
                label faceStart = pf.patch().start();

                forAll(pf, facei)
                {
                    allValues[faceStart++] = pf[facei];
                }
            }

            this->distribute(allValues);
            newValues.transfer(allValues);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown sampling mode: " << this->mode() << nl
                << abort(FatalError);
        }
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_ && returnReduceOr(newValues.size()))
    {
        Type averagePsi;
        if (this->faceValues())
        {
            const scalarField magSf
            (
                mag(this->mappedPatchBase::patch_.faceAreas())
            );
            averagePsi = gSum(magSf*newValues)/gSum(magSf);
        }
        else
        {
            averagePsi = gAverage(newValues);
        }

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

    return this->transform(tnewValues);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::Sampled<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
void Foam::PatchFunction1Types::Sampled<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);

    os.writeEntry(this->name(), type());

    mappedPatchBase::write(os);

    os.writeEntry("field", fieldName_);
    if (setAverage_)
    {
        os.writeEntry("setAverage", "true");
        os.writeEntry("average", average_);
    }

    os.writeEntry("interpolationScheme", interpolationScheme_);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "outletMappedUniformInletFvPatchField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValuePtr_(nullptr),
    outletNames_(),
    offsets_(),
    fractions_(),
    timeDelays_(),
    mapFields_(),
    mapTimes_(),
    phiName_("phi"),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    uniformValuePtr_
    (
        PatchFunction1<Type>::NewIfPresent
        (
            p.patch(),
            "uniformValue",
            dict
        )
    ),
    outletNames_(),
    offsets_(),
    fractions_(),
    timeDelays_(),
    mapFields_(),
    mapTimes_(),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    curTimeIndex_(-1)
{
    const dictionary& outletDict = dict.subDict("outlets");

    if (outletDict.empty())
    {
        FatalIOErrorInFunction(outletDict)
            << "outlets dictionary is empty."
            << exit(FatalIOError);
    }

    outletNames_.setSize(outletDict.size());
    offsets_.setSize(outletDict.size());
    fractions_.setSize(outletDict.size());
    timeDelays_.setSize(outletDict.size());
    mapFields_.setSize(outletDict.size());
    mapTimes_.setSize(outletDict.size());

    label outleti = 0;
    for (const entry& dEntry : outletDict)
    {
        const word& key = dEntry.keyword();

        if (!dEntry.isDict())
        {
            FatalIOErrorInFunction(outletDict)
                << "Entry " << key << " is not a dictionary." << nl
                << exit(FatalIOError);
        }

        const dictionary& subDict = dEntry.dict();

        outletNames_[outleti] = key;

        offsets_.set
        (
            outleti,
            Function1<Type>::NewIfPresent
            (
                "offset",
                subDict,
                word::null,
                &this->db()
            )
        );

        fractions_.set
        (
            outleti,
            Function1<scalar>::NewIfPresent
            (
                "fraction",
                subDict,
                word::null,
                &this->db()
            )
        );

        timeDelays_.set
        (
            outleti,
            Function1<scalar>::NewIfPresent
            (
                "timeDelay",
                subDict,
                word::null,
                &this->db()
            )
        );

        mapFields_[outleti] =
            subDict.getOrDefault<DynamicList<Type>>
            (
                "mapField",
                DynamicList<Type>()
            );

        mapTimes_[outleti] =
            subDict.getOrDefault<DynamicList<scalar>>
            (
                "mapTime",
                DynamicList<scalar>()
            );

        ++outleti;
    }


    if (!this->readValueEntry(dict))
    {
        // Fallback: set to the internal field
        this->extrapolateInternal();
    }
}


template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const outletMappedUniformInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    uniformValuePtr_(ptf.uniformValuePtr_.clone(p.patch())),
    outletNames_(ptf.outletNames_),
    offsets_(ptf.offsets_),
    fractions_(ptf.fractions_),
    timeDelays_(ptf.timeDelays_),
    mapFields_(ptf.mapFields_),
    mapTimes_(ptf.mapTimes_),
    phiName_(ptf.phiName_),
    curTimeIndex_(-1)
{
    if (mapper.direct() && !mapper.hasUnmapped())
    {
        // Use mapping instead of re-evaluation
        this->map(ptf, mapper);
    }
    else
    {
        // Fallback: set to the internal field
        this->extrapolateInternal();
    }
}


template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const outletMappedUniformInletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValuePtr_(ptf.uniformValuePtr_.clone(this->patch().patch())),
    outletNames_(ptf.outletNames_),
    offsets_(ptf.offsets_),
    fractions_(ptf.fractions_),
    timeDelays_(ptf.timeDelays_),
    mapFields_(ptf.mapFields_),
    mapTimes_(ptf.mapTimes_),
    phiName_(ptf.phiName_),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::outletMappedUniformInletFvPatchField<Type>::
outletMappedUniformInletFvPatchField
(
    const outletMappedUniformInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValuePtr_(ptf.uniformValuePtr_.clone(this->patch().patch())),
    outletNames_(ptf.outletNames_),
    offsets_(ptf.offsets_),
    fractions_(ptf.fractions_),
    timeDelays_(ptf.timeDelays_),
    mapFields_(ptf.mapFields_),
    mapTimes_(ptf.mapTimes_),
    phiName_(ptf.phiName_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::outletMappedUniformInletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);

    if (uniformValuePtr_)
    {
        uniformValuePtr_->autoMap(m);
    }
}


template<class Type>
void Foam::outletMappedUniformInletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const auto& tiptf =
        refCast<const outletMappedUniformInletFvPatchField>(ptf);

    if (uniformValuePtr_)
    {
        uniformValuePtr_->rmap(tiptf.uniformValuePtr_(), addr);
    }
}


template<class Type>
void Foam::outletMappedUniformInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        const scalar t = this->db().time().timeOutputValue();

        const GeometricField<Type, fvPatchField, volMesh>& f
        (
            dynamic_cast<const GeometricField<Type, fvPatchField, volMesh>&>
            (
                this->internalField()
            )
        );

        const fvPatch& p = this->patch();

        forAll(outletNames_, i)
        {
            const word& outletName = outletNames_[i];
            const label outletID =
                p.patch().boundaryMesh().findPatchID(outletName);

            if (outletID < 0)
            {
                FatalErrorInFunction
                    << "Unable to find outlet patch " << outletName
                    << abort(FatalError);
            }


            // Collect the map time for this outlet patch
            DynamicList<scalar>& mapTime = mapTimes_[i];
            scalar timeDelay = 0;
            if (timeDelays_.set(i))
            {
                timeDelay = max(timeDelays_[i].value(t), scalar(0));
            }
            mapTime.append(t + timeDelay);


            // Collect the map field for this outlet patch and map time
            const fvPatchField<Type>& outletFld = f.boundaryField()[outletID];
            DynamicList<Type>& mapField = mapFields_[i];

            const auto& phi =
                this->db().objectRegistry::template
                lookupObject<surfaceScalarField>(phiName_);
            const scalarField& outletPhi = phi.boundaryField()[outletID];
            const scalar sumOutletPhi = gSum(outletPhi);

            if (sumOutletPhi > SMALL)
            {
                Type offset(Zero);
                if (offsets_.set(i))
                {
                    offset = offsets_[i].value(t);
                }

                scalar fraction = 1;
                if (fractions_.set(i))
                {
                    fraction = fractions_[i].value(t);
                }

                mapField.append
                (
                    gSum(outletPhi*outletFld)/sumOutletPhi*fraction
                  + offset
                );
            }
            else
            {
                const fvPatch& outlet = p.boundaryMesh()[outletID];

                mapField.append
                (
                    gSum(outlet.magSf()*outletFld)/gSum(outlet.magSf())
                );
            }
        }


        // Map the stored fields onto inlet if the time condition is met
        Type inletFld(Zero);
        forAll(outletNames_, i)
        {
            DynamicList<scalar>& mapTime = mapTimes_[i];
            DynamicList<Type>& mapField = mapFields_[i];

            if (!mapTime.empty())
            {
                if (t >= mapTime.first())
                {
                    inletFld += interpolateXY(t, mapTime, mapField);

                    // Remove any stored fields and times if possible
                    int i = 0;
                    while (!mapTime.empty() && t >= mapTime[i])
                    {
                        mapTime.remove(i);
                        mapField.remove(i);
                        ++i;
                    }
                }
            }
        }


        if (uniformValuePtr_)
        {
            this->operator==(inletFld + uniformValuePtr_->value(t));
        }
        else
        {
            this->operator==(inletFld);
        }


        curTimeIndex_ = this->db().time().timeIndex();
    }


    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::outletMappedUniformInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    if (uniformValuePtr_)
    {
        uniformValuePtr_->writeData(os);
    }
    os.beginBlock("outlets");
    forAll(outletNames_, i)
    {
        os.beginBlock(outletNames_[i]);
        if (offsets_.set(i))
        {
            offsets_[i].writeData(os);
        }
        if (fractions_.set(i))
        {
            fractions_[i].writeData(os);
        }
        if (timeDelays_.set(i))
        {
            timeDelays_[i].writeData(os);
        }
        if (!mapFields_.empty())
        {
            mapFields_[i].writeEntry("mapField", os);
        }
        if (!mapTimes_.empty())
        {
            mapTimes_[i].writeEntry("mapTime", os);
        }
        os.endBlock();
    }
    os.endBlock();
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);

    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "timeVaryingMappedFixedValuePointPatchField.H"
#include "Time.H"
#include "rawIOField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMappedFixedValuePointPatchField<Type>::
timeVaryingMappedFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    setAverage_(false),
    perturb_(0),
    fieldTableName_(iF.name()),
    pointsName_("points"),
    mapMethod_(),
    mapperPtr_(nullptr),
    sampleTimes_(),
    sampleIndex_(-1, -1),
    sampleAverage_(Zero, Zero),
    sampleValues_(),
    offset_(nullptr)
{}


template<class Type>
Foam::timeVaryingMappedFixedValuePointPatchField<Type>::
timeVaryingMappedFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    setAverage_(dict.getOrDefault("setAverage", false)),
    perturb_(dict.getOrDefault("perturb", 1e-5)),
    fieldTableName_(iF.name()),
    pointsName_(dict.getOrDefault<word>("points", "points")),
    mapMethod_(),
    mapperPtr_(nullptr),
    sampleTimes_(),
    sampleIndex_(-1, -1),
    sampleAverage_(Zero, Zero),
    sampleValues_(),
    offset_
    (
        Function1<Type>::NewIfPresent("offset", dict, word::null, &this->db())
    )
{
    if
    (
        dict.readIfPresent("mapMethod", mapMethod_)
     && !mapMethod_.empty()
     && mapMethod_ != "nearest"
     && !mapMethod_.starts_with("planar")
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown mapMethod type " << mapMethod_
            << "\n\nValid mapMethod types :\n"
            << "(nearest planar)" << nl
            << exit(FatalIOError);
    }

    dict.readIfPresentCompat
    (
        "fieldTable", {{"fieldTableName", 2206}},
        fieldTableName_
    );

    if (!this->readValueEntry(dict))
    {
        // Note: use evaluate to do updateCoeffs followed by a reset
        //       of the pointPatchField::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        pointPatchField<Type>::evaluate(Pstream::commsTypes::buffered);
    }
}


template<class Type>
Foam::timeVaryingMappedFixedValuePointPatchField<Type>::
timeVaryingMappedFixedValuePointPatchField
(
    const timeVaryingMappedFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    fieldTableName_(ptf.fieldTableName_),
    pointsName_(ptf.pointsName_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(nullptr),
    sampleTimes_(),
    sampleIndex_(-1, -1),
    sampleAverage_(Zero, Zero),
    sampleValues_(),
    offset_(ptf.offset_.clone())
{}


template<class Type>
Foam::timeVaryingMappedFixedValuePointPatchField<Type>::
timeVaryingMappedFixedValuePointPatchField
(
    const timeVaryingMappedFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    fieldTableName_(ptf.fieldTableName_),
    pointsName_(ptf.pointsName_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(ptf.mapperPtr_),
    sampleTimes_(ptf.sampleTimes_),
    sampleIndex_(ptf.sampleIndex_),
    sampleAverage_(ptf.sampleAverage_),
    sampleValues_(ptf.sampleValues_),
    offset_(ptf.offset_.clone())
{}


template<class Type>
Foam::timeVaryingMappedFixedValuePointPatchField<Type>::
timeVaryingMappedFixedValuePointPatchField
(
    const timeVaryingMappedFixedValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    setAverage_(ptf.setAverage_),
    perturb_(ptf.perturb_),
    fieldTableName_(ptf.fieldTableName_),
    pointsName_(ptf.pointsName_),
    mapMethod_(ptf.mapMethod_),
    mapperPtr_(ptf.mapperPtr_),
    sampleTimes_(ptf.sampleTimes_),
    sampleIndex_(ptf.sampleIndex_),
    sampleAverage_(ptf.sampleAverage_),
    sampleValues_(ptf.sampleValues_),
    offset_(ptf.offset_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedFixedValuePointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    fixedValuePointPatchField<Type>::autoMap(m);

    if (sampleValues_.first().size())
    {
        sampleValues_.first().autoMap(m);
    }

    if (sampleValues_.second().size())
    {
        sampleValues_.second().autoMap(m);
    }

    // Clear interpolator
    mapperPtr_.reset(nullptr);
    sampleIndex_ = labelPair(-1, -1);
}


template<class Type>
void Foam::timeVaryingMappedFixedValuePointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValuePointPatchField<Type>::rmap(ptf, addr);

    const timeVaryingMappedFixedValuePointPatchField<Type>& tiptf =
        refCast<const timeVaryingMappedFixedValuePointPatchField<Type>>(ptf);

    sampleValues_.first().rmap(tiptf.sampleValues_.first(), addr);
    sampleValues_.second().rmap(tiptf.sampleValues_.second(), addr);

    // Clear interpolator
    mapperPtr_.reset(nullptr);
    sampleIndex_ = labelPair(-1, -1);
}


template<class Type>
void Foam::timeVaryingMappedFixedValuePointPatchField<Type>::updateSampledValues
(
    const label sampleIndex,
    Field<Type>& field,
    Type& avg
) const
{
    tmp<Field<Type>> tvalues;

    // Update sampled data fields
    {
        const Time& time = this->db().time();

        if (debug)
        {
            Pout<< "checkTable : Reading values from "
                <<
                (
                    "boundaryData"
                  / this->patch().name()
                  / sampleTimes_[sampleIndex].name()
                  / fieldTableName_
                ) << endl;
        }

        // Reread values and interpolate
        const fileName valsFile
        (
            time.caseConstant()
            /"boundaryData"
            /this->patch().name()
            /sampleTimes_[sampleIndex].name()
            /fieldTableName_
        );

        IOobject io
        (
            valsFile,               // absolute path
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER,
            true                    // is global object (currently not used)
        );

        rawIOField<Type> vals(io, setAverage_);

        if (vals.hasAverage())
        {
            avg = vals.average();
        }

        if (vals.size() != mapperPtr_().sourceSize())
        {
            FatalErrorInFunction
                << "Number of values (" << vals.size()
                << ") differs from the number of points ("
                <<  mapperPtr_().sourceSize()
                << ") in file " << valsFile << exit(FatalError);
        }

        tvalues = tmp<Field<Type>>::New(std::move(vals.field()));
    }

    // From input values to interpolated (sampled) positions
    field = mapperPtr_().interpolate(tvalues);
}


template<class Type>
void Foam::timeVaryingMappedFixedValuePointPatchField<Type>::checkTable
(
    const scalar t
)
{
   const Time& time = this->db().time();

    // Initialise
    if (sampleIndex_.first() == -1 && sampleIndex_.second() == -1)
    {
        const polyMesh& pMesh = this->patch().boundaryMesh().mesh()();

        // Read the initial point position
        pointField meshPts;

        if (pMesh.pointsInstance() == pMesh.facesInstance())
        {
            meshPts = pointField(pMesh.points(), this->patch().meshPoints());
        }
        else
        {
            // Load points from facesInstance
            DebugInFunction
                << "Reloading points0 from " << pMesh.facesInstance()
                << endl;

            pointIOField points0
            (
                IOobject
                (
                    "points",
                    pMesh.facesInstance(),
                    polyMesh::meshSubDir,
                    pMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );
            meshPts = pointField(points0, this->patch().meshPoints());
        }

        // Reread values and interpolate
        const fileName samplePointsFile
        (
            time.caseConstant()
           /"boundaryData"
           /this->patch().name()
           /pointsName_
        );

        IOobject io
        (
            samplePointsFile,   // absolute path
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER,
            true                // is global object (currently not used)
        );

        // Read data (no average value!)
        const rawIOField<point> samplePoints(io);

        // tbd: run-time selection
        const bool nearestOnly =
        (
            !mapMethod_.empty() && !mapMethod_.starts_with("planar")
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                meshPts,
                perturb_,
                nearestOnly
            )
        );


        // Read the times for which data is available
        const fileName samplePointsDir = samplePointsFile.path();
        sampleTimes_ = Time::findTimes(samplePointsDir);

        DebugInFunction
            << "Found times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_)
            << endl;
    }

    // Find range of current time indices in sampleTimes
    Pair<label> timeIndices = instant::findRange
    (
        sampleTimes_,
        t,  // time.value(),
        sampleIndex_.first()
    );

    if (timeIndices.first() < 0)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for current time "
            << time.value() << nl
            << "Have sampling values for times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  time.constant()/"boundaryData"/this->patch().name()
            << "\n    on patch " << this->patch().name()
            << " of field " << fieldTableName_
            << exit(FatalError);
    }


    // Update sampled data fields.

    if (sampleIndex_.first() != timeIndices.first())
    {
        sampleIndex_.first() = timeIndices.first();

        if (sampleIndex_.first() == sampleIndex_.second())
        {
            // Can reuse previous end values
            sampleValues_.first() = sampleValues_.second();
            sampleAverage_.first() = sampleAverage_.second();
        }
        else
        {
            // Update first() values
            this->updateSampledValues
            (
                sampleIndex_.first(),
                sampleValues_.first(),
                sampleAverage_.first()
            );
        }
    }

    if (sampleIndex_.second() != timeIndices.second())
    {
        sampleIndex_.second() = timeIndices.second();

        if (sampleIndex_.second() == -1)
        {
            // Index no longer valid - clear values accordingly
            sampleValues_.second().clear();
        }
        else
        {
            // Update second() values
            this->updateSampledValues
            (
                sampleIndex_.second(),
                sampleValues_.second(),
                sampleAverage_.second()
            );
        }
    }
}


template<class Type>
void Foam::timeVaryingMappedFixedValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Current time value
    const scalar x = this->db().time().value();

    checkTable(x);

    // Interpolate between the sampled data
    auto& fld = static_cast<Field<Type>&>(*this);
    Type wantedAverage;

    if (sampleIndex_.second() == -1)
    {
        // Only start value
        fld = sampleValues_.first();
        wantedAverage = sampleAverage_.first();
    }
    else
    {
        const scalar beg = sampleTimes_[sampleIndex_.first()].value();
        const scalar end = sampleTimes_[sampleIndex_.second()].value();
        const scalar s = (x - beg)/(end - beg);

        fld = (1 - s)*sampleValues_.first() + s*sampleValues_.second();

        wantedAverage
            = (1 - s)*sampleAverage_.first() + s*sampleAverage_.second();

        DebugInFunction
            << "Sampled, interpolated values"
            << " between time:"
            << sampleTimes_[sampleIndex_.first()].name()
            << " and time:" << sampleTimes_[sampleIndex_.second()].name()
            << " with weight:" << s << endl;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        Type averagePsi = gAverage(fld);

        if (debug)
        {
            Pout<< "updateCoeffs :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < VSMALL)
        {
            // Field too small to scale. Offset instead.
            const Type offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " offsetting with:" << offset << endl;
            }
            fld += offset;
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " scaling with:" << scale << endl;
            }
            fld *= scale;
        }
    }

    // Apply offset to mapped values
    if (offset_)
    {
        const scalar t = this->db().time().timeOutputValue();
        fld += offset_->value(t);
    }

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this) << endl;
    }

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingMappedFixedValuePointPatchField<Type>::write
(
    Ostream& os
) const
{
    fixedValuePointPatchField<Type>::write(os);

    os.writeEntryIfDifferent
    (
        "fieldTable",
        this->internalField().name(),
        fieldTableName_
    );

    if (!pointsName_.empty())
    {
        os.writeEntryIfDifferent<word>("points", "points", pointsName_);
    }

    if (!mapMethod_.empty() && !mapMethod_.starts_with("planar"))
    {
        os.writeEntry("mapMethod", mapMethod_);
    }

    if (setAverage_)
    {
        os.writeEntry("setAverage", setAverage_);
    }

    os.writeEntryIfDifferent<scalar>("perturb", 1e-5, perturb_);

    if (offset_)
    {
        offset_->writeData(os);
    }
}


// ************************************************************************* //

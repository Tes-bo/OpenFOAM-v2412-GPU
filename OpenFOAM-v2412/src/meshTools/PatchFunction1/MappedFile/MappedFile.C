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

#include "polyMesh.H"
#include "rawIOField.H"
#include "clockTime.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const bool dictConstructed,
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const word& fieldTableName,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    dictConstructed_(dictConstructed),
    setAverage_(dict.getOrDefault("setAverage", false)),
    perturb_(dict.getOrDefault<scalar>("perturb", 1e-5)),
    fieldTableName_(fieldTableName),
    pointsName_(dict.getOrDefault<word>("points", "points")),
    mapMethod_(),
    filterRadius_(dict.getOrDefault<scalar>("filterRadius", 0)),
    filterSweeps_(dict.getOrDefault<label>("filterSweeps", 0)),
    filterFieldPtr_(nullptr),
    readerFormat_(),
    readerFile_(),
    readerPtr_(nullptr),
    mapperPtr_(nullptr),
    sampleTimes_(),
    sampleIndex_(-1, -1),
    sampleAverage_(Zero, Zero),
    sampleValues_(),
    offset_
    (
        Function1<Type>::NewIfPresent
        (
            "offset",
            dict,
            patchFunction1Base::whichDb()
        )
    )
{
    if (fieldTableName_.empty())
    {
        fieldTableName_ = entryName;
    }

    // Simple sanity check
    if ((filterSweeps_ < 1) || (filterRadius_ <= VSMALL))
    {
        filterRadius_ = 0;
        filterSweeps_ = 0;
    }

    if (dict.readIfPresent("sampleFormat", readerFormat_))
    {
        dict.readEntry("sampleFile", readerFile_);

        fileName fName(readerFile_);
        fName.expand();

        readerPtr_ = surfaceReader::New
        (
            readerFormat_,
            fName,
            surfaceReader::formatOptions(dict, readerFormat_, "readOptions")
        );
    }

    if (debug)
    {
        Info<< "mappedFile:" << nl;
        if (readerFormat_.empty())
        {
            Info<< "    boundary format" << nl;
        }
        else
        {
            Info<< "    format:" << readerFormat_
                << " file:" << readerFile_ << nl;
        }

        Info<< "    filter radius=" << filterRadius_
            << " sweeps=" << filterSweeps_ << endl;
    }

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
}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    MappedFile<Type>
    (
        true,  // dictConstructed = true
        pp,
        entryName,
        dict,
        dict.getOrDefault<word>("fieldTable", entryName),
        faceValues
    )
{}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const word& fieldTableName,
    const bool faceValues
)
:
    MappedFile<Type>
    (
        false,  // dictConstructed = false
        pp,
        entryName,
        dict,
        fieldTableName,
        faceValues
    )
{}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const MappedFile<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    dictConstructed_(rhs.dictConstructed_),
    setAverage_(rhs.setAverage_),
    perturb_(rhs.perturb_),
    fieldTableName_(rhs.fieldTableName_),
    pointsName_(rhs.pointsName_),
    mapMethod_(rhs.mapMethod_),
    filterRadius_(rhs.filterRadius_),
    filterSweeps_(rhs.filterSweeps_),
    filterFieldPtr_(nullptr),
    readerFormat_(rhs.readerFormat_),
    readerFile_(rhs.readerFile_),
    readerPtr_(nullptr),
    mapperPtr_(rhs.mapperPtr_.clone()),
    sampleTimes_(rhs.sampleTimes_),
    sampleIndex_(rhs.sampleIndex_),
    sampleAverage_(rhs.sampleAverage_),
    sampleValues_(rhs.sampleValues_),
    offset_(rhs.offset_.clone())
{
    if (!readerFormat_.empty() && !readerFile_.empty())
    {
        fileName fName(readerFile_);
        fName.expand();

        readerPtr_ = surfaceReader::New(readerFormat_, fName);
    }
}


template<class Type>
Foam::PatchFunction1Types::MappedFile<Type>::MappedFile
(
    const MappedFile<Type>& rhs
)
:
    MappedFile<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::autoMap
(
    const FieldMapper& mapper
)
{
    PatchFunction1<Type>::autoMap(mapper);

    if (sampleValues_.first().size())
    {
        sampleValues_.first().autoMap(mapper);
    }
    if (sampleValues_.second().size())
    {
        sampleValues_.second().autoMap(mapper);
    }

    // Clear interpolator
    filterFieldPtr_.reset(nullptr);
    mapperPtr_.reset(nullptr);
    sampleIndex_ = labelPair(-1, -1);
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    PatchFunction1<Type>::rmap(pf1, addr);

    const auto& tiptf =
        refCast<const PatchFunction1Types::MappedFile<Type>>(pf1);

    if (tiptf.sampleValues_.first().size())
    {
        sampleValues_.first().resize(this->size());
        sampleValues_.first().rmap(tiptf.sampleValues_.first(), addr);
    }

    if (tiptf.sampleValues_.second().size())
    {
        sampleValues_.second().resize(this->size());
        sampleValues_.second().rmap(tiptf.sampleValues_.second(), addr);
    }

    // Clear interpolator
    filterFieldPtr_.reset(nullptr);
    mapperPtr_.reset(nullptr);
    sampleIndex_ = labelPair(-1, -1);
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::updateSampledValues
(
    const label sampleIndex,
    Field<Type>& field,
    Type& avg
) const
{
    tmp<Field<Type>> tvalues;

    // Update sampled data fields
    if (readerPtr_)
    {
        wordList fieldNames = readerPtr_->fieldNames(sampleIndex);

        label fieldIndex = fieldNames.find(fieldTableName_);

        if (fieldIndex < 0)
        {
            FatalErrorInFunction
                << "Sample field='" << fieldTableName_
                << "' not found. Known field names: "
                << flatOutput(fieldNames) << nl
                << exit(FatalError);
        }

        if (debug)
        {
            Pout<< "checkTable : Update index=" << sampleIndex
                << " field=" << fieldNames[fieldIndex] << endl;
        }

        tvalues = readerPtr_->field
        (
            sampleIndex,
            fieldIndex,
            pTraits<Type>::zero
        );

        if (tvalues().size() != mapperPtr_().sourceSize())
        {
            FatalErrorInFunction
                << "Number of values (" << tvalues().size()
                << ") differs from the number of points ("
                <<  mapperPtr_().sourceSize() << ")"
                << exit(FatalError);
        }
    }
    else
    {
        const polyMesh& mesh = this->patch_.boundaryMesh().mesh();
        const Time& time = mesh.time();

        if (debug)
        {
            Pout<< "checkTable : Update index=" << sampleIndex
                << " Reading values from "
                <<
                (
                    "boundaryData"
                  / this->patch_.name()
                  / sampleTimes_[sampleIndex].name()
                  / fieldTableName_
                ) << endl;
        }

        // Reread values and interpolate
        const fileName valsFile
        (
            time.globalPath()
            /time.constant()
            /mesh.dbDir()            // region
            /"boundaryData"
            /this->patch_.name()
            /sampleTimes_[sampleIndex].name()
            /fieldTableName_
        );

        IOobject io
        (
            valsFile,           // absolute path
            time,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER,
            true                // is global object (currently not used)
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
                << ") in file " << valsFile
                << exit(FatalError);
        }

        tvalues = tmp<Field<Type>>::New(std::move(vals.field()));
    }

    if (filterFieldPtr_)
    {
        DebugInfo
            << "apply " << filterSweeps_ << " filter sweeps" << endl;

        tvalues = filterFieldPtr_().evaluate(tvalues, filterSweeps_);
    }

    // From input values to interpolated (sampled) positions
    field = mapperPtr_().interpolate(tvalues);
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::checkTable
(
    const scalar t
) const
{
    const polyMesh& mesh = this->patch_.boundaryMesh().mesh();
    const Time& time = mesh.time();

    // Initialise
    if (!mapperPtr_ && readerPtr_)
    {
        clockTime timing;

        auto& reader = readerPtr_();

        const meshedSurface& geom = reader.geometry(0);

        sampleTimes_ = reader.times();

        // We may have been passed in geometry without any faces
        // eg, boundaryData

        const pointField& samplePoints =
        (
            geom.nFaces() ? geom.faceCentres() : geom.points()
        );

        DebugInfo
            << "Read " << samplePoints.size() << " sample points from "
            << readerFile_ << endl
            << "Found times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "... in " << timing.timeIncrement() << 's' << endl;

        // tbd: run-time selection
        const bool nearestOnly =
        (
            !mapMethod_.empty() && !mapMethod_.starts_with("planar")
        );

        // Allocate the interpolator
        if (this->faceValues())
        {
            mapperPtr_.reset
            (
                new pointToPointPlanarInterpolation
                (
                    samplePoints,
                    this->localPosition(this->patch_.faceCentres()),
                    perturb_,
                    nearestOnly
                )
            );
        }
        else
        {
            mapperPtr_.reset
            (
                new pointToPointPlanarInterpolation
                (
                    samplePoints,
                    this->localPosition(this->patch_.localPoints()),
                    perturb_,
                    nearestOnly
                )
            );
        }

        DebugInfo
            << "Created point/point planar interpolation"
            << " - in " << timing.timeIncrement() << 's' << endl;


        // Setup median filter (if any)
        if (filterSweeps_ > 0)
        {
            filterFieldPtr_.reset(new FilterField(geom, filterRadius_));

            DebugInfo
                << "Calculated field-filter"
                << " - in " << timing.timeIncrement() << 's' << endl;
        }
        else
        {
            filterFieldPtr_.reset(nullptr);
        }
    }
    else if (!mapperPtr_)
    {
        clockTime timing;

        // Reread values and interpolate
        const fileName samplePointsFile
        (
            time.globalPath()
           /time.constant()         // instance
           /mesh.dbDir()            // region
           /"boundaryData"
           /this->patch_.name()
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

        // Read the times for which data is available
        sampleTimes_ = Time::findTimes(samplePointsFile.path());

        DebugInfo
            << "Read " << samplePoints.size() << " sample points from "
            << samplePointsFile << endl
            << "Found times "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "... in " << timing.timeIncrement() << 's' << endl;


        // tbd: run-time selection
        const bool nearestOnly =
        (
            !mapMethod_.empty() && !mapMethod_.starts_with("planar")
        );

        // Allocate the interpolator
        if (this->faceValues())
        {
            mapperPtr_.reset
            (
                new pointToPointPlanarInterpolation
                (
                    samplePoints,
                    this->localPosition(this->patch_.faceCentres()),
                    perturb_,
                    nearestOnly
                )
            );
        }
        else
        {
            mapperPtr_.reset
            (
                new pointToPointPlanarInterpolation
                (
                    samplePoints,
                    this->localPosition(this->patch_.localPoints()),
                    perturb_,
                    nearestOnly
                )
            );
        }

        DebugInfo
            << "Created point/point planar interpolation"
            << " - in " << timing.timeIncrement() << 's' << endl;


        // Setup median filter (if any)
        if (filterSweeps_ > 0)
        {
            filterFieldPtr_.reset(new FilterField(samplePoints, filterRadius_));

            DebugInfo
                << "Calculated field-filter"
                << " - in " << timing.timeIncrement() << 's' << endl;
        }
        else
        {
            filterFieldPtr_.reset(nullptr);
        }
    }


    // Find range of current time indices in sampleTimes
    labelPair timeIndices = instant::findRange
    (
        sampleTimes_,
        t,  //mesh.time().value(),
        sampleIndex_.first()
    );

    if (timeIndices.first() < 0)
    {
        FatalErrorInFunction
            << "Cannot find starting sampling values for index "
            << t << nl
            << "Have sampling values for "
            << pointToPointPlanarInterpolation::timeNames(sampleTimes_) << nl
            << "In directory "
            <<  time.constant()/mesh.dbDir()/"boundaryData"/this->patch_.name()
            << "\n    on patch " << this->patch_.name()
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
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::MappedFile<Type>::value
(
    const scalar x
) const
{
    checkTable(x);

    // Interpolate between the sampled data
    auto tfld = tmp<Field<Type>>::New();
    auto& fld = tfld.ref();
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

        DebugInfo
            << "MappedFile<Type>::value : "
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
        Type averagePsi;

        if (this->faceValues())
        {
            const scalarField magSf(mag(this->patch_.faceAreas()));
            averagePsi = gSum(magSf*fld)/gSum(magSf);
        }
        else
        {
            averagePsi = gAverage(fld);
        }

        if (debug)
        {
            Pout<< "MappedFile<Type>::value :"
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
                Pout<< "MappedFile<Type>::value :"
                    << " offsetting with:" << offset << endl;
            }
            fld += offset;
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "MappedFile<Type>::value :"
                    << " scaling with:" << scale << endl;
            }
            fld *= scale;
        }
    }

    // Apply offset to mapped values
    if (offset_)
    {
        fld += offset_->value(x);
    }

    if (debug)
    {
        Pout<< "MappedFile<Type>::value : set fixedValue to min:" << gMin(fld)
            << " max:" << gMax(fld)
            << " avg:" << gAverage(fld) << endl;
    }

    return this->transform(tfld);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::MappedFile<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::writeEntries
(
    Ostream& os
) const
{
    if (!readerFormat_.empty() && !readerFile_.empty())
    {
        os.writeEntry("readerFormat", readerFormat_);
        os.writeEntry("readerFile", readerFile_);
    }

    os.writeEntryIfDifferent
    (
        "fieldTable",
        this->name(),
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

    if (filterSweeps_ >= 1)
    {
        os.writeEntry("filterRadius", filterRadius_);
        os.writeEntry("filterSweeps", filterSweeps_);
    }

    if (offset_)
    {
        offset_->writeData(os);
    }
}


template<class Type>
void Foam::PatchFunction1Types::MappedFile<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);

    // Check if field name explicitly provided
    // (e.g. through timeVaryingMapped bc)
    if (dictConstructed_)
    {
        os.writeEntry(this->name(), type());

        os.beginBlock(word(this->name() + "Coeffs"));
        writeEntries(os);
        os.endBlock();
    }
    else
    {
        // Note that usually dictConstructed = true. The
        // construct-from-dictionary (dictConstructed_ = false)
        // should only be used if there is only
        // a single potential MappedFile in the local scope.
        writeEntries(os);
    }
}


// ************************************************************************* //

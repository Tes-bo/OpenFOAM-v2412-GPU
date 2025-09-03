/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "surfaceWriter.H"
#include "proxySurfaceWriter.H"
#include "MeshedSurfaceProxy.H"

#include "fileFormats.H"
#include "Time.H"
#include "coordinateRotation.H"
#include "transformField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceWriter, 0);
    defineRunTimeSelectionTable(surfaceWriter, word);
    defineRunTimeSelectionTable(surfaceWriter, wordDict);
}

Foam::scalar Foam::surfaceWriter::defaultMergeDim = 1e-8;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::surfaceWriter::supportedType(const word& writeType)
{
    return
    (
        wordConstructorTablePtr_->found(writeType)
     || wordDictConstructorTablePtr_->found(writeType)
     || MeshedSurfaceProxy<face>::canWriteType(writeType)
    );
}


Foam::dictionary Foam::surfaceWriter::formatOptions
(
    const dictionary& dict,
    const word& formatName,
    const word& entryName
)
{
    return fileFormats::getFormatOptions(dict, formatName, entryName);
}


Foam::dictionary Foam::surfaceWriter::formatOptions
(
    const dictionary& dict,
    const dictionary& surfDict,
    const word& formatName,
    const word& entryName
)
{
    return fileFormats::getFormatOptions(dict, surfDict, formatName, entryName);
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::TryNew(const word& writeType)
{
    // Without dictionary options
    {
        auto* ctorPtr = wordConstructorTable(writeType);

        if (ctorPtr)
        {
            return autoPtr<surfaceWriter>(ctorPtr());
        }
    }

    // Fallback to proxy writer...
    return surfaceWriters::proxyWriter::TryNew(writeType);
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::TryNew
(
    const word& writeType,
    const dictionary& writeOpts
)
{
    // With dictionary options
    {
        auto* ctorPtr = wordDictConstructorTable(writeType);

        if (ctorPtr)
        {
            return autoPtr<surfaceWriter>(ctorPtr(writeOpts));
        }
    }

    // Without dictionary options
    {
        auto* ctorPtr = wordConstructorTable(writeType);

        if (ctorPtr)
        {
            return autoPtr<surfaceWriter>(ctorPtr());
        }
    }

    // Fallback to proxy writer...
    return surfaceWriters::proxyWriter::TryNew(writeType, writeOpts);
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New(const word& writeType)
{
    autoPtr<surfaceWriter> writer
    (
        surfaceWriter::TryNew(writeType)
    );

    if (!writer)
    {
        FatalErrorInFunction
            << "Unknown write type \"" << writeType << "\"\n\n"
            << "Valid write types : "
            << flatOutput(wordConstructorTablePtr_->sortedToc()) << nl
            << "Valid proxy types : "
            << flatOutput(MeshedSurfaceProxy<face>::writeTypes()) << endl
            << exit(FatalError);
    }

    return writer;
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New
(
    const word& writeType,
    const dictionary& writeOpts
)
{
    autoPtr<surfaceWriter> writer
    (
        surfaceWriter::TryNew(writeType, writeOpts)
    );

    if (!writer)
    {
        FatalErrorInFunction
            << "Unknown write type \"" << writeType << "\"\n\n"
            << "Valid write types : "
            << flatOutput(wordConstructorTablePtr_->sortedToc()) << nl
            << "Valid proxy types : "
            << flatOutput(MeshedSurfaceProxy<face>::writeTypes()) << endl
            << exit(FatalError);
    }

    return writer;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriter::surfaceWriter()
:
    surf_(),
    mergedSurf_(),
    adjustedSurf_(),
    mergeDim_(defaultMergeDim),
    geometryScale_(1),
    geometryCentre_(Zero),
    geometryTransform_(),
    upToDate_(false),
    wroteGeom_(false),
    parallel_(true),
    useTimeDir_(false),
    isPointData_(false),
    verbose_(false),
    commType_(UPstream::commsTypes::scheduled),
    nFields_(0),
    currTime_(),
    outputPath_(),
    fieldLevel_(),
    fieldScale_()
{
    surfaceWriter::close();
}


Foam::surfaceWriter::surfaceWriter(const dictionary& options)
:
    surfaceWriter()
{
    options.readIfPresent("verbose", verbose_);

    UPstream::commsTypeNames.readIfPresent("commsType", options, commType_);

    geometryScale_ = 1;
    geometryCentre_ = Zero;
    geometryTransform_.clear();

    options.readIfPresent("scale", geometryScale_);

    // Optional cartesian coordinate system transform
    const auto* dictptr = options.findDict("transform", keyType::LITERAL);

    if (dictptr)
    {
        dictptr->readIfPresent("rotationCentre", geometryCentre_);

        // 'origin' is optional within sub-dictionary
        geometryTransform_ =
            coordSystem::cartesian(*dictptr, IOobjectOption::LAZY_READ);
    }

    fieldLevel_ = options.subOrEmptyDict("fieldLevel");
    fieldScale_ = options.subOrEmptyDict("fieldScale");

    if (verbose_)
    {
        Info<< "Create surfaceWriter ("
            << (this->isPointData() ? "point" : "face") << " data):"
            << " commsType=" << UPstream::commsTypeNames[commType_] << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceWriter::~surfaceWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceWriter::setTime(const instant& inst)
{
    currTime_ = inst;
}


void Foam::surfaceWriter::setTime(scalar timeValue)
{
    currTime_ = instant(timeValue);
}


void Foam::surfaceWriter::setTime(scalar timeValue, const word& timeName)
{
    currTime_.value() = timeValue;
    currTime_.name() = timeName;
}


void Foam::surfaceWriter::unsetTime()
{
    currTime_.value() = 0;
    currTime_.name().clear();
}


void Foam::surfaceWriter::beginTime(const Time& t)
{
    setTime(t.value(), t.timeName());
}


void Foam::surfaceWriter::beginTime(const instant& inst)
{
    setTime(inst);
}


void Foam::surfaceWriter::endTime()
{
    unsetTime();
}


void Foam::surfaceWriter::open(const fileName& outputPath)
{
    outputPath_ = outputPath;
    wroteGeom_ = false;
}


void Foam::surfaceWriter::open
(
    const meshedSurf& surf,
    const fileName& outputPath,
    bool parallel
)
{
    close();
    setSurface(surf, parallel);
    open(outputPath);
}


void Foam::surfaceWriter::open
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath,
    bool parallel
)
{
    close();
    setSurface(points, faces, parallel);
    open(outputPath);
}


void Foam::surfaceWriter::open
(
    const meshedSurf& surf,
    const fileName& outputPath
)
{
    close();
    setSurface(surf, parallel_);
    open(outputPath);
}


void Foam::surfaceWriter::open
(
    const pointField& points,
    const faceList& faces,
    const fileName& outputPath
)
{
    close();
    setSurface(points, faces, parallel_);
    open(outputPath);
}


void Foam::surfaceWriter::close()
{
    outputPath_.clear();
    wroteGeom_ = false;
}


void Foam::surfaceWriter::clear()
{
    close();
    expire();
    surf_.clear();
}


void Foam::surfaceWriter::setSurface
(
    const meshedSurf& surf,
    bool parallel
)
{
    expire();
    surf_.reset(surf);
    parallel_ = (parallel && UPstream::parRun());
}


void Foam::surfaceWriter::setSurface
(
    const pointField& points,
    const faceList& faces,
    bool parallel
)
{
    expire();
    surf_.reset(points, faces);
    parallel_ = (parallel && UPstream::parRun());
}


void Foam::surfaceWriter::setSurface
(
    const meshedSurf& surf
)
{
    setSurface(surf, parallel_);
}


void Foam::surfaceWriter::setSurface
(
    const pointField& points,
    const faceList& faces
)
{
    setSurface(points, faces, parallel_);
}


bool Foam::surfaceWriter::needsUpdate() const
{
    return !upToDate_;
}


bool Foam::surfaceWriter::wroteData() const
{
    return wroteGeom_;
}


bool Foam::surfaceWriter::expire()
{
    const bool changed = upToDate_;

    upToDate_ = false;
    wroteGeom_ = false;
    adjustedSurf_.clear();
    mergedSurf_.clear();

    // Field count (nFields_) is a different type of accounting
    // and is unaffected by geometry changes

    return changed;
}


bool Foam::surfaceWriter::hasSurface() const
{
    return surf_.good();
}


bool Foam::surfaceWriter::empty() const
{
    const bool value = surf_.faces().empty();

    return (parallel_ ? returnReduceAnd(value) : value);
}


Foam::label Foam::surfaceWriter::size() const
{
    const label value = surf_.faces().size();

    return (parallel_ ? returnReduce(value, sumOp<label>()) : value);
}


void Foam::surfaceWriter::checkOpen() const
{
    if (!is_open())
    {
        FatalErrorInFunction
            << type() << " : Attempted to write without a path" << nl
            << exit(FatalError);
    }
}


bool Foam::surfaceWriter::merge() const
{
    bool changed = false;

    if (!upToDate_)
    {
        // Similar to expire
        adjustedSurf_.clear();

        if (parallel_ && UPstream::parRun())
        {
            changed = mergedSurf_.merge(surf_, mergeDim_);
        }
        else
        {
            mergedSurf_.clear();
        }

        if (changed)
        {
            wroteGeom_ = false;
        }
    }

    upToDate_ = true;
    return changed;
}


const Foam::meshedSurf& Foam::surfaceWriter::surface() const
{
    merge();

    if (parallel_ && UPstream::parRun())
    {
        return mergedSurf_;
    }

    return surf_;
}


const Foam::meshedSurfRef& Foam::surfaceWriter::adjustSurface() const
{
    if (!upToDate_)
    {
        adjustedSurf_.clear();
    }

    if (!adjustedSurf_.good())
    {
        adjustedSurf_.reset(surface());

        tmp<pointField> tpts;

        if (geometryTransform_.good())
        {
            if (!geometryTransform_.R().is_identity())
            {
                if (magSqr(geometryCentre_) > ROOTVSMALL)
                {
                    // Set centre of rotation,
                    // followed by forward transform (local -> global)
                    tpts =
                        geometryTransform_.globalPosition
                        (
                            adjustedSurf_.points0() - geometryCentre_
                        );

                    // Unset centre of rotation
                    tpts.ref() += geometryCentre_;
                }
                else
                {
                    // Forward transform (local -> global)
                    tpts =
                        geometryTransform_.globalPosition
                        (
                            adjustedSurf_.points0()
                        );
                }
            }
            else if (magSqr(geometryTransform_.origin()) > ROOTVSMALL)
            {
                // Translate only (local -> global)
                tpts = (adjustedSurf_.points0() + geometryTransform_.origin());
            }
        }

        adjustedSurf_.movePoints(tpts);
        adjustedSurf_.scalePoints(geometryScale_);
    }

    return adjustedSurf_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::surfaceWriter::mergeFieldTemplate
(
    const Field<Type>& fld
) const
{
    if (parallel_ && UPstream::parRun())
    {
        // Ensure geometry is also merged
        merge();

        // Gather all values
        auto tfield = tmp<Field<Type>>::New();
        auto& allFld = tfield.ref();

        // Update any expired global index (as required)

        const globalIndex& globIndex =
        (
            this->isPointData()
          ? mergedSurf_.pointGlobalIndex()
          : mergedSurf_.faceGlobalIndex()
        );

        globIndex.gather
        (
            fld,
            allFld,
            UPstream::msgType(),
            commType_,
            UPstream::worldComm
        );

        // Renumber (point data) to correspond to merged points
        if
        (
            UPstream::master()
         && this->isPointData()
         && mergedSurf_.pointsMap().size()
        )
        {
            inplaceReorder(mergedSurf_.pointsMap(), allFld);
            allFld.resize(mergedSurf_.points().size());
        }

        return tfield;
    }

    // Mark that any geometry changes have been taken care of
    upToDate_ = true;

    return fld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::surfaceWriter::adjustFieldTemplate
(
    const word& fieldName,
    const tmp<Field<Type>>& tfield
) const
{
    if (verbose_)
    {
        Info<< "Writing field " << fieldName;
    }

    tmp<Field<Type>> tadjusted;

    // Output scaling for the variable, but not for integer types
    // which are typically ids etc.
    if (!std::is_integral<Type>::value)
    {
        scalar value;

        // Remove *uniform* reference level
        if
        (
            fieldLevel_.readIfPresent(fieldName, value, keyType::REGEX)
         && !equal(value, 0)
        )
        {
            // Could also detect brackets (...) and read accordingly
            // or automatically scale by 1/sqrt(nComponents) instead ...

            Type refLevel;
            for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
            {
                setComponent(refLevel, cmpt) = value;
            }

            if (verbose_)
            {
                Info<< " [level " << refLevel << ']';
            }

            if (!tadjusted)
            {
                // Steal or clone
                tadjusted.reset(tfield.ptr());
            }

            // Remove offset level
            tadjusted.ref() -= refLevel;
        }

        // Apply scaling
        if
        (
            fieldScale_.readIfPresent(fieldName, value, keyType::REGEX)
         && !equal(value, 1)
        )
        {
            if (verbose_)
            {
                Info<< " [scaling " << value << ']';
            }

            if (!tadjusted)
            {
                // Steal or clone
                tadjusted.reset(tfield.ptr());
            }

            // Apply scaling
            tadjusted.ref() *= value;
        }

        // Rotate fields (vector and non-spherical tensors)
        if
        (
            (is_vectorspace<Type>::value && pTraits<Type>::nComponents > 1)
         && geometryTransform_.good()
         && !geometryTransform_.R().is_identity()
        )
        {
            if (!tadjusted)
            {
                // Steal or clone
                tadjusted.reset(tfield.ptr());
            }

            Foam::transform
            (
                tadjusted.ref(),
                geometryTransform_.R(),
                tadjusted()
            );
        }
    }

    return (tadjusted ? tadjusted : tfield);
}


#define defineSurfaceFieldMethods(ThisClass, Type)                             \
    Foam::tmp<Foam::Field<Type>>                                               \
    ThisClass::mergeField(const Field<Type>& fld) const                        \
    {                                                                          \
        return mergeFieldTemplate(fld);                                        \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::Field<Type>>                                               \
    ThisClass::adjustField                                                     \
    (                                                                          \
        const word& fieldName,                                                 \
        const tmp<Field<Type>>& tfield                                         \
    ) const                                                                    \
    {                                                                          \
        return adjustFieldTemplate(fieldName, tfield);                         \
    }

defineSurfaceFieldMethods(Foam::surfaceWriter, Foam::label);
defineSurfaceFieldMethods(Foam::surfaceWriter, Foam::scalar);
defineSurfaceFieldMethods(Foam::surfaceWriter, Foam::vector);
defineSurfaceFieldMethods(Foam::surfaceWriter, Foam::sphericalTensor);
defineSurfaceFieldMethods(Foam::surfaceWriter, Foam::symmTensor);
defineSurfaceFieldMethods(Foam::surfaceWriter, Foam::tensor)

#undef defineSurfaceFieldMethod


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<surfaceWriter>& iproxy
)
{
    const auto& w = *iproxy;

    os  << "surfaceWriter:"
        << " upToDate: " << w.upToDate_
        << " PointData: " << w.isPointData_
        << " nFields: " << w.nFields_
        << " time: " << w.currTime_
        << " path: " << w.outputPath_ << endl;

    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "surfaceFieldValue.H"
#include "fvMesh.H"
#include "emptyPolyPatch.H"
#include "coupledPolyPatch.H"
#include "sampledSurface.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Max number of warnings
static constexpr const unsigned maxWarnings = 10u;

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(surfaceFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, surfaceFieldValue, runTime);
    addToRunTimeSelectionTable(functionObject, surfaceFieldValue, dictionary);
}
}
}


const Foam::Enum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypes
>
Foam::functionObjects::fieldValues::surfaceFieldValue::regionTypeNames_
({
    { regionTypes::stFaceZone, "faceZone" },
    { regionTypes::stPatch, "patch" },
    { regionTypes::stObject, "functionObjectSurface" },
    { regionTypes::stSampled, "sampledSurface" },
});


const Foam::Enum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::operationType
>
Foam::functionObjects::fieldValues::surfaceFieldValue::operationTypeNames_
({
    // Normal operations
    { operationType::opNone, "none" },
    { operationType::opMin, "min" },
    { operationType::opMax, "max" },
    { operationType::opSum, "sum" },
    { operationType::opSumMag, "sumMag" },
    { operationType::opSumDirection, "sumDirection" },
    { operationType::opSumDirectionBalance, "sumDirectionBalance" },
    { operationType::opAverage, "average" },
    { operationType::opAreaAverage, "areaAverage" },
    { operationType::opAreaIntegrate, "areaIntegrate" },
    { operationType::opCoV, "CoV" },
    { operationType::opAreaNormalAverage, "areaNormalAverage" },
    { operationType::opAreaNormalIntegrate, "areaNormalIntegrate" },
    { operationType::opUniformity, "uniformity" },

    // Using weighting
    { operationType::opWeightedSum, "weightedSum" },
    { operationType::opWeightedAverage, "weightedAverage" },
    { operationType::opWeightedAreaAverage, "weightedAreaAverage" },
    { operationType::opWeightedAreaIntegrate, "weightedAreaIntegrate" },
    { operationType::opWeightedUniformity, "weightedUniformity" },

    // Using absolute weighting
    { operationType::opAbsWeightedSum, "absWeightedSum" },
    { operationType::opAbsWeightedAverage, "absWeightedAverage" },
    { operationType::opAbsWeightedAreaAverage, "absWeightedAreaAverage" },
    { operationType::opAbsWeightedAreaIntegrate, "absWeightedAreaIntegrate" },
    { operationType::opAbsWeightedUniformity, "absWeightedUniformity" },
});

const Foam::Enum
<
    Foam::functionObjects::fieldValues::surfaceFieldValue::postOperationType
>
Foam::functionObjects::fieldValues::surfaceFieldValue::postOperationTypeNames_
({
    { postOperationType::postOpNone, "none" },
    { postOperationType::postOpMag, "mag" },
    { postOperationType::postOpSqrt, "sqrt" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::objectRegistry&
Foam::functionObjects::fieldValues::surfaceFieldValue::obr() const
{
    if (stObject == regionType_)
    {
        return storedObjects().lookupObject<polySurface>(regionName_);
    }

    return mesh_;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setFaceZoneFaces()
{
    // Indices for all matches, already sorted
    const labelList zoneIds
    (
        mesh_.faceZones().indices(selectionNames_)
    );

    // Total number of faces that could be selected (before patch filtering)
    label numFaces = 0;
    for (const label zoneId : zoneIds)
    {
        numFaces += mesh_.faceZones()[zoneId].size();
    }

    faceId_.resize_nocopy(numFaces);
    facePatchId_.resize_nocopy(numFaces);
    faceFlip_.resize_nocopy(numFaces);

    numFaces = 0;

    for (const label zoneId : zoneIds)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneId];

        forAll(fZone, i)
        {
            const label meshFacei = fZone[i];
            const bool isFlip = fZone.flipMap()[i];

            // Internal faces
            label faceId = meshFacei;
            label facePatchId = -1;

            // Boundary faces
            if (!mesh_.isInternalFace(meshFacei))
            {
                facePatchId = mesh_.boundaryMesh().whichPatch(meshFacei);
                const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];

                if (isA<emptyPolyPatch>(pp))
                {
                    continue;  // Ignore empty patch
                }

                const auto* cpp = isA<coupledPolyPatch>(pp);

                if (cpp && !cpp->owner())
                {
                    continue;  // Ignore neighbour side
                }

                faceId = pp.whichFace(meshFacei);
            }

            if (faceId >= 0)
            {
                faceId_[numFaces] = faceId;
                facePatchId_[numFaces] = facePatchId;
                faceFlip_[numFaces] = isFlip;

                ++numFaces;
            }
        }
    }

    // Shrink to size used
    faceId_.resize(numFaces);
    facePatchId_.resize(numFaces);
    faceFlip_.resize(numFaces);
    nFaces_ = returnReduce(numFaces, sumOp<label>());

    if (!nFaces_)
    {
        // Raise warning or error
        refPtr<OSstream> os;
        bool fatal = false;

        ++nWarnings_;  // Always increment (even if ignore etc)

        switch (emptySurfaceError_)
        {
            case error::handlerTypes::IGNORE:
            {
                break;
            }

            case error::handlerTypes::WARN:
            {
                if (nWarnings_ <= maxWarnings)
                {
                    os.ref(WarningInFunction);
                }
                break;
            }

            // STRICT / DEFAULT
            default:
            {
                os.ref(FatalErrorInFunction);
                fatal = true;
                break;
            }
        }

        if (os)
        {
            os.ref()
                << type() << ' ' << name() << ": "
                << regionTypeNames_[regionType_]
                << '(' << regionName_ << "):" << nl;

            if (zoneIds.empty())
            {
                os.ref()
                    << "    No matching face zones: "
                    << flatOutput(selectionNames_)  << nl
                    << "    Known face zones: "
                    << flatOutput(mesh_.faceZones().names()) << nl;
            }
            else
            {
                os.ref()
                    << "    The face zones: "
                    << flatOutput(selectionNames_) << nl
                    << "    resulted in 0 faces" << nl;
            }

            if (fatal)
            {
                FatalError<< exit(FatalError);
            }
            else if (nWarnings_ == maxWarnings)
            {
                os.ref()
                    << "... suppressing further warnings." << nl;
            }
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::setPatchFaces()
{
    // Patch indices for all matches
    labelList patchIds;

    // Total number of faces selected
    label numFaces = 0;

    labelList selected
    (
        mesh_.boundaryMesh().indices(selectionNames_, true)  // useGroup
    );

    DynamicList<label> bad;
    for (const label patchi : selected)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        if (isA<emptyPolyPatch>(pp))
        {
            bad.append(patchi);
        }
        else
        {
            numFaces += pp.size();
        }
    }

    if (bad.size())
    {
        label nGood = (selected.size() - bad.size());

        auto& os = (nGood > 0 ? WarningInFunction : FatalErrorInFunction);

        os  << "Cannot sample an empty patch" << nl;

        for (const label patchi : bad)
        {
            os  << "    "
                << mesh_.boundaryMesh()[patchi].name() << nl;
        }

        if (nGood)
        {
            os  << "No non-empty patches selected" << endl
                << exit(FatalError);
        }
        else
        {
            os  << "Selected " << nGood << " non-empty patches" << nl;
        }

        patchIds.resize(nGood);
        nGood = 0;
        for (const label patchi : selected)
        {
            if (!bad.contains(patchi))
            {
                patchIds[nGood] = patchi;
                ++nGood;
            }
        }
    }
    else
    {
        patchIds = std::move(selected);
    }

    faceId_.resize_nocopy(numFaces);
    facePatchId_.resize_nocopy(numFaces);
    faceFlip_.resize_nocopy(numFaces);
    nFaces_ = returnReduce(numFaces, sumOp<label>());

    numFaces = 0;
    for (const label patchi : patchIds)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        const label len = pp.size();

        SubList<label>(faceId_, len, numFaces) = identity(len);
        SubList<label>(facePatchId_, len, numFaces) = patchi;
        SubList<bool>(faceFlip_, len, numFaces) = false;

        numFaces += len;
    }

    if (!nFaces_)
    {
        // Raise warning or error
        refPtr<OSstream> os;
        bool fatal = false;

        ++nWarnings_;  // Always increment (even if ignore etc)

        switch (emptySurfaceError_)
        {
            case error::handlerTypes::IGNORE:
            {
                break;
            }

            case error::handlerTypes::WARN:
            {
                if (nWarnings_ <= maxWarnings)
                {
                    os.ref(WarningInFunction);
                }
                break;
            }

            // STRICT / DEFAULT
            default:
            {
                os.ref(FatalErrorInFunction);
                fatal = true;
                break;
            }
        }

        if (os)
        {
            os.ref()
                << type() << ' ' << name() << ": "
                << regionTypeNames_[regionType_]
                << '(' << regionName_ << "):" << nl;

            if (patchIds.empty())
            {
                os.ref()
                    << "    No matching patches: "
                    << flatOutput(selectionNames_)  << nl
                    << "    Known patch names:" << nl
                    << mesh_.boundaryMesh().names() << nl;
            }
            else
            {
                os.ref()
                    << "    The patches: "
                    << flatOutput(selectionNames_) << nl
                    << "    resulted in 0 faces" << nl;
            }

            if (fatal)
            {
                FatalError<< exit(FatalError);
            }
            else if (nWarnings_ == maxWarnings)
            {
                os.ref()
                    << "... suppressing further warnings." << nl;
            }
        }
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::combineMeshGeometry
(
    faceList& faces,
    pointField& points
) const
{
    labelList whichFaces(faceId_);

    // Remap patch-face ids to mesh face ids
    forAll(whichFaces, i)
    {
        const label patchi = facePatchId_[i];
        if (patchi != -1)
        {
            whichFaces[i] += mesh_.boundaryMesh()[patchi].start();
        }
    }

    indirectPrimitivePatch pp
    (
        IndirectList<face>(mesh_.faces(), std::move(whichFaces)),
        mesh_.points()
    );


    if (Pstream::parRun())
    {
        // Topological merge
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPoints;
        autoPtr<globalIndex> globalFaces;
        Foam::PatchTools::gatherAndMerge
        (
            mesh_,
            pp.localFaces(),
            pp.meshPoints(),
            pp.meshPointMap(),

            pointToGlobal,
            uniqueMeshPointLabels,
            globalPoints,
            globalFaces,

            faces,
            points
        );
    }
    else
    {
        faces = pp.localFaces();
        points = pp.localPoints();
    }
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::
combineSurfaceGeometry
(
    faceList& faces,
    pointField& points
) const
{
    if (stObject == regionType_)
    {
        const auto& s = refCast<const polySurface>(obr());

        if (Pstream::parRun())
        {
            // Dimension as fraction of surface
            const scalar mergeDim = 1e-10*boundBox(s.points(), true).mag();

            Foam::PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch(SubList<face>(s.faces()), s.points()),
                points,
                faces
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
    else if (sampledPtr_)
    {
        const sampledSurface& s = *sampledPtr_;

        if (Pstream::parRun())
        {
            // Dimension as fraction of mesh bounding box
            const scalar mergeDim = 1e-10*mesh_.bounds().mag();

            Foam::PatchTools::gatherAndMerge
            (
                mergeDim,
                primitivePatch(SubList<face>(s.faces()), s.points()),
                points,
                faces
            );
        }
        else
        {
            faces = s.faces();
            points = s.points();
        }
    }
}


Foam::scalar
Foam::functionObjects::fieldValues::surfaceFieldValue::totalArea() const
{
    scalar totalArea = 0;

    if (stObject == regionType_)
    {
        const auto& s = refCast<const polySurface>(obr());

        totalArea = gSum(s.magSf());
    }
    else if (sampledPtr_)
    {
        totalArea = gSum(sampledPtr_->magSf());
    }
    else
    {
        totalArea = gSum(filterField(mesh_.magSf()));
    }

    return totalArea;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValue::usesSf()
const noexcept
{
    // Few operations do not require the Sf field
    switch (operation_)
    {
        case opNone:
        case opMin:
        case opMax:
        case opSum:
        case opSumMag:
        case opAverage:
            return false;

        default:
            return true;
    }
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::update()
{
    if (sampledPtr_)
    {
        sampledPtr_->update();
    }

    if (!needsUpdate_)
    {
        return false;
    }

    // Reset some values
    totalArea_ = 0;
    nFaces_ = 0;
    bool checkEmptyFaces = true;

    switch (regionType_)
    {
        case stFaceZone:
        {
            // Raises warning or error internally, don't check again
            setFaceZoneFaces();
            checkEmptyFaces = false;
            break;
        }
        case stPatch:
        {
            // Raises warning or error internally, don't check again
            setPatchFaces();
            checkEmptyFaces = false;
            break;
        }
        case stObject:
        {
            // TBD: special handling of cast errors?
            const auto& s = refCast<const polySurface>(obr());
            nFaces_ = returnReduce(s.size(), sumOp<label>());
            break;
        }
        case stSampled:
        {
            nFaces_ = returnReduce(sampledPtr_->faces().size(), sumOp<label>());
            break;
        }

        // Compiler warning if we forgot an enumeration
    }

    if (nFaces_)
    {
        // Appears to be successful
        needsUpdate_ = false;
        totalArea_ = totalArea();   // Update the area
        nWarnings_ = 0u;            // Reset the warnings counter
    }
    else if (checkEmptyFaces)
    {
        // Raise warning or error
        refPtr<OSstream> os;
        bool fatal = false;

        ++nWarnings_;  // Always increment (even if ignore etc)

        switch (emptySurfaceError_)
        {
            case error::handlerTypes::IGNORE:
            {
                break;
            }

            case error::handlerTypes::WARN:
            {
                if (nWarnings_ <= maxWarnings)
                {
                    os.ref(WarningInFunction);
                }
                break;
            }

            // STRICT / DEFAULT
            default:
            {
                os.ref(FatalErrorInFunction);
                fatal = true;
                break;
            }
        }

        if (os)
        {
            os.ref()
                << type() << ' ' << name() << ": "
                << regionTypeNames_[regionType_]
                << '(' << regionName_ << "):" << nl
                << "    Region has no faces" << endl;

            if (fatal)
            {
                FatalError<< exit(FatalError);
            }
            else if (nWarnings_ == maxWarnings)
            {
                os.ref()
                    << "... suppressing further warnings." << nl;
            }
        }
    }

    Log << "    total faces   = " << nFaces_ << nl
        << "    total area    = " << totalArea_ << nl
        << endl;

    // Emit file header on success or change of state
    if (nWarnings_ <= 1)
    {
        writeFileHeader(file());
    }

    return true;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::writeFileHeader
(
    Ostream& os
)
{
    if (canWriteHeader() && (operation_ != opNone))
    {
        writeCommented(os, "Region type : ");
        os << regionTypeNames_[regionType_] << ' ' << regionName_ << nl;

        writeHeaderValue(os, "Faces", nFaces_);
        writeHeaderValue(os, "Area", totalArea_);
        writeHeaderValue(os, "Scale factor", scaleFactor_);

        if (weightFieldNames_.size())
        {
            writeHeaderValue
            (
                os,
                "Weight field",
                flatOutput(weightFieldNames_, FlatOutput::BareComma{})
            );
        }

        writeCommented(os, "Time");
        if (writeArea_)
        {
            os  << tab << "Area";
        }

        // TBD: add in postOperation information?

        for (const word& fieldName : fields_)
        {
            os  << tab << operationTypeNames_[operation_]
                << '(' << fieldName << ')';
        }

        os  << endl;
    }

    writtenHeader_ = true;
}


template<>
Foam::scalar
Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<scalar>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            const vector n(dict_.get<vector>("direction"));
            return gSum(pos0(values*(Sf & n))*mag(values));
        }
        case opSumDirectionBalance:
        {
            const vector n(dict_.get<vector>("direction"));
            const scalarField nv(values*(Sf & n));

            return gSum(pos0(nv)*mag(values) - neg(nv)*mag(values));
        }

        case opUniformity:
        case opWeightedUniformity:
        case opAbsWeightedUniformity:
        {
            const scalar areaTotal = gSum(mag(Sf));
            tmp<scalarField> areaVal(values * mag(Sf));

            scalar mean, numer;

            if (is_weightedOp() && canWeight(weightField))
            {
                // Weighted quantity = (Weight * phi * dA)

                tmp<scalarField> weight
                (
                    weightingFactor(weightField, is_magOp())
                );

                // Mean weighted value (area-averaged)
                mean = gSum(weight()*areaVal()) / areaTotal;

                // Abs. deviation from weighted mean value
                numer = gSum(mag(weight*areaVal - (mean * mag(Sf))));
            }
            else
            {
                // Unweighted quantity = (1 * phi * dA)

                // Mean value (area-averaged)
                mean = gSum(areaVal()) / areaTotal;

                // Abs. deviation from mean value
                numer = gSum(mag(areaVal - (mean * mag(Sf))));
            }

            // Uniformity index
            const scalar ui = 1 - numer/(2*mag(mean*areaTotal) + ROOTVSMALL);

            return clamp(ui, zero_one{});
        }

        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


template<>
Foam::vector
Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<vector>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const
{
    switch (operation_)
    {
        case opSumDirection:
        {
            const vector n(dict_.get<vector>("direction").normalise());

            const scalarField nv(n & values);
            return gSum(pos0(nv)*n*(nv));
        }
        case opSumDirectionBalance:
        {
            const vector n(dict_.get<vector>("direction").normalise());

            const scalarField nv(n & values);
            return gSum(pos0(nv)*n*(nv));
        }
        case opAreaNormalAverage:
        {
            const scalar val = gSum(values & Sf)/gSum(mag(Sf));
            return vector(val, 0, 0);
        }
        case opAreaNormalIntegrate:
        {
            const scalar val = gSum(values & Sf);
            return vector(val, 0, 0);
        }

        case opUniformity:
        case opWeightedUniformity:
        case opAbsWeightedUniformity:
        {
            const scalar areaTotal = gSum(mag(Sf));
            tmp<scalarField> areaVal(values & Sf);

            scalar mean, numer;

            if (is_weightedOp() && canWeight(weightField))
            {
                // Weighted quantity = (Weight * phi . dA)

                tmp<scalarField> weight
                (
                    weightingFactor(weightField, is_magOp())
                );

                // Mean weighted value (area-averaged)
                mean = gSum(weight()*areaVal()) / areaTotal;

                // Abs. deviation from weighted mean value
                numer = gSum(mag(weight*areaVal - (mean * mag(Sf))));
            }
            else
            {
                // Unweighted quantity = (1 * phi . dA)

                // Mean value (area-averaged)
                mean = gSum(areaVal()) / areaTotal;

                // Abs. deviation from mean value
                numer = gSum(mag(areaVal - (mean * mag(Sf))));
            }

            // Uniformity index
            const scalar ui = 1 - numer/(2*mag(mean*areaTotal) + ROOTVSMALL);

            return vector(clamp(ui, zero_one{}), 0, 0);
        }

        default:
        {
            // Fall through to other operations
            return processSameTypeValues(values, Sf, weightField);
        }
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField,
    const bool useMag
)
{
    if (useMag)
    {
        return mag(weightField);
    }

    // pass through
    return weightField;
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField,
    const vectorField& Sf,
    const bool useMag
)
{
    // scalar * unit-normal

    // Can skip this check - already used canWeight()
    /// if (returnReduceAnd(weightField.empty()))
    /// {
    ///     // No weight field - revert to unweighted form?
    ///     return tmp<scalarField>::New(Sf.size(), scalar(1));
    /// }

    if (useMag)
    {
        return mag(weightField);
    }

    // pass through
    return weightField;
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::areaWeightingFactor
(
    const Field<scalar>& weightField,
    const vectorField& Sf,
    const bool useMag
)
{
    // scalar * Area

    // Can skip this check - already used canWeight()
    /// if (returnReduceAnd(weightField.empty()))
    /// {
    ///     // No weight field - revert to unweighted form
    ///     return mag(Sf);
    /// }

    if (useMag)
    {
        return mag(weightField * mag(Sf));
    }

    return (weightField * mag(Sf));
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<vector>& weightField,
    const vectorField& Sf,
    const bool useMag
)
{
    // vector (dot) unit-normal

    // Can skip this check - already used canWeight()
    /// if (returnReduceAnd(weightField.empty()))
    /// {
    ///     // No weight field - revert to unweighted form
    ///     return tmp<scalarField>::New(Sf.size(), scalar(1));
    /// }

    const label len = weightField.size();

    auto tresult = tmp<scalarField>::New(weightField.size());
    auto& result = tresult.ref();

    for (label facei=0; facei < len; ++facei)
    {
        const vector unitNormal(normalised(Sf[facei]));
        result[facei] = (weightField[facei] & unitNormal);
    }

    if (useMag)
    {
        for (scalar& val : result)
        {
            val = mag(val);
        }
    }

    return tresult;
}


template<>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::areaWeightingFactor
(
    const Field<vector>& weightField,
    const vectorField& Sf,
    const bool useMag
)
{
    // vector (dot) Area

    // Can skip this check - already used canWeight()
    /// if (returnReduceAnd(weightField.empty()))
    /// {
    ///     // No weight field - revert to unweighted form
    ///     return mag(Sf);
    /// }

    if (useMag)
    {
        return mag(weightField & Sf);
    }

    return (weightField & Sf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    regionType_(regionTypeNames_.get("regionType", dict)),
    operation_(operationTypeNames_.get("operation", dict)),
    postOperation_
    (
        postOperationTypeNames_.getOrDefault
        (
            "postOperation",
            dict,
            postOperationType::postOpNone,
            true  // Failsafe behaviour
        )
    ),
    needsUpdate_(true),
    writeArea_(false),
    emptySurfaceError_(error::handlerTypes::DEFAULT),
    selectionNames_(),
    weightFieldNames_(),
    totalArea_(0),
    nFaces_(0),
    nWarnings_(0),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);
}


Foam::functionObjects::fieldValues::surfaceFieldValue::surfaceFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    regionType_(regionTypeNames_.get("regionType", dict)),
    operation_(operationTypeNames_.get("operation", dict)),
    postOperation_
    (
        postOperationTypeNames_.getOrDefault
        (
            "postOperation",
            dict,
            postOperationType::postOpNone,
            true  // Failsafe behaviour
        )
    ),
    needsUpdate_(true),
    writeArea_(false),
    emptySurfaceError_(error::handlerTypes::DEFAULT),
    selectionNames_(),
    weightFieldNames_(),
    totalArea_(0),
    nFaces_(0),
    nWarnings_(0),
    faceId_(),
    facePatchId_(),
    faceFlip_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Needs completed sampledSurface, surfaceWriter
Foam::functionObjects::fieldValues::surfaceFieldValue::~surfaceFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::surfaceFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);

    needsUpdate_ = true;
    writeArea_ = dict.getOrDefault("writeArea", false);
    emptySurfaceError_ = error::handlerNames.getOrDefault
    (
        "empty-surface",
        dict,
        error::handlerTypes::DEFAULT,
        true  // Failsafe behaviour
    );

    weightFieldNames_.clear();
    // future?
    // sampleFaceScheme_ = dict.getOrDefault<word>("sampleScheme", "cell");

    totalArea_ = 0;
    nFaces_ = 0;
    nWarnings_ = 0;
    faceId_.clear();
    facePatchId_.clear();
    faceFlip_.clear();
    sampledPtr_.reset(nullptr);
    surfaceWriterPtr_.reset(nullptr);

    // Can have "name" (word) and/or "names" (wordRes)
    //
    // If "names" exists AND contains a literal (non-regex) that can be used
    // as a suitable value for "name", the "name" entry becomes optional.

    regionName_.clear();
    selectionNames_.clear();

    {
        dict.readIfPresent("names", selectionNames_);

        for (const auto& item : selectionNames_)
        {
            if (item.isLiteral())
            {
                regionName_ = item;
                break;
            }
        }

        // The "name" entry
        // - mandatory if we didn't pick up a value from selectionNames_
        dict.readEntry
        (
            "name",
            regionName_,
            keyType::LITERAL,
            (
                regionName_.empty()
              ? IOobjectOption::MUST_READ : IOobjectOption::READ_IF_PRESENT
            )
        );

        // Ensure there is always content for selectionNames_
        if (selectionNames_.empty())
        {
            selectionNames_.resize(1);
            selectionNames_.first() = regionName_;
        }
    }


    // Create sampled surface, but leave 'expired' (ie, no update) since it
    // may depend on fields or data that do not yet exist
    if (stSampled == regionType_)
    {
        sampledPtr_ = sampledSurface::New
        (
            name(),
            mesh_,
            dict.subDict("sampledSurfaceDict")
        );

        // Internal consistency. Want face values, never point values!
        sampledPtr_->isPointData(false);
    }

    Info<< type() << ' ' << name() << ':' << nl
        << "    operation     = ";

    if (postOperation_ != postOpNone)
    {
        Info<< postOperationTypeNames_[postOperation_] << '('
            << operationTypeNames_[operation_] << ')'  << nl;
    }
    else
    {
        Info<< operationTypeNames_[operation_] << nl;
    }

    if (is_weightedOp())
    {
        // Can have "weightFields" or "weightField"

        bool missing = true;
        if (dict.readIfPresent("weightFields", weightFieldNames_))
        {
            missing = false;
        }
        else
        {
            weightFieldNames_.resize(1);

            if (dict.readIfPresent("weightField", weightFieldNames_.first()))
            {
                missing = false;
                if ("none" == weightFieldNames_.first())
                {
                    // "none" == no weighting
                    weightFieldNames_.clear();
                }
            }
        }

        if (missing)
        {
            // Suggest possible alternative unweighted operation?
            FatalIOErrorInFunction(dict)
                << "The '" << operationTypeNames_[operation_]
                << "' operation is missing a weightField." << nl
                << "Either provide the weightField, "
                << "use weightField 'none' to suppress weighting," << nl
                << "or use a different operation."
                << exit(FatalIOError);
        }

        Info<< "    weight field  = ";
        if (weightFieldNames_.empty())
        {
            Info<< "none" << nl;
        }
        else
        {
            Info<< flatOutput(weightFieldNames_) << nl;
        }
    }

    if (stSampled == regionType_ && sampledPtr_)
    {
        Info<< "    sampled surface: ";
        sampledPtr_->print(Info, 0);
        Info<< nl;
    }

    if (writeFields_)
    {
        const word writerType = dict.get<word>("surfaceFormat");

        surfaceWriterPtr_.reset
        (
            surfaceWriter::New
            (
                writerType,
                surfaceWriter::formatOptions(dict, writerType)
            )
        );

        // Propagate field counts (per surface)
        surfaceWriterPtr_->nFields(fields_.size());

        if (debug)
        {
            surfaceWriterPtr_->verbose(true);
        }

        if (surfaceWriterPtr_->enabled())
        {
            Info<< "    surfaceFormat = " << writerType << nl;
        }
        else
        {
            surfaceWriterPtr_->clear();
        }
    }

    Info<< nl << endl;

    return true;
}


bool Foam::functionObjects::fieldValues::surfaceFieldValue::write()
{
    if (needsUpdate_ || operation_ != opNone)
    {
        fieldValue::write();
    }

    update();

    if (operation_ != opNone)
    {
        writeCurrentTime(file());
    }

    // Handle ignore/warn about empty-surface
    if (!nFaces_)
    {
        totalArea_ = 0;  // Update the area (safety)

        if (operation_ != opNone)
        {
            if (emptySurfaceError_ == error::handlerTypes::WARN)
            {
                if (writeArea_)
                {
                    Log << "    total area = " << totalArea_ << endl;
                    file() << tab << totalArea_;
                }

                file() << tab << "NaN";
                Log << endl;
            }

            file() << endl;
        }

        // Early exit on error
        return true;
    }

    if (writeArea_)
    {
        // Update the area
        totalArea_ = totalArea();
        Log << "    total area = " << totalArea_ << endl;

        if (operation_ != opNone && UPstream::master())
        {
            file() << tab << totalArea_;
        }
    }

    // Many operations use the Sf field
    vectorField Sf;
    if (usesSf())
    {
        if (stObject == regionType_)
        {
            const auto& s = refCast<const polySurface>(obr());
            Sf = s.Sf();
        }
        else if (sampledPtr_)
        {
            Sf = sampledPtr_->Sf();
        }
        else
        {
            Sf = filterField(mesh_.Sf());
        }
    }

    // Faces and points for surface format (if specified)
    faceList faces;
    pointField points;

    if (surfaceWriterPtr_)
    {
        if (withTopologicalMerge())
        {
            combineMeshGeometry(faces, points);
        }
        else
        {
            combineSurfaceGeometry(faces, points);
        }
    }


    // Check availability and type of weight field
    // Only support a few weight types:
    // scalar: 0-N fields
    // vector: 0-1 fields

    // Default is a zero-size scalar weight field (ie, weight = 1)
    scalarField scalarWeights;
    vectorField vectorWeights;

    for (const word& weightName : weightFieldNames_)
    {
        if (validField<scalar>(weightName))
        {
            tmp<scalarField> tfld = getFieldValues<scalar>(weightName, true);

            if (scalarWeights.empty())
            {
                scalarWeights = tfld;
            }
            else
            {
                scalarWeights *= tfld;
            }
        }
        else if (validField<vector>(weightName))
        {
            tmp<vectorField> tfld = getFieldValues<vector>(weightName, true);

            if (vectorWeights.empty())
            {
                vectorWeights = tfld;
            }
            else
            {
                FatalErrorInFunction
                    << "weightField " << weightName
                    << " - only one vector weight field allowed. " << nl
                    << "weights: " << flatOutput(weightFieldNames_) << nl
                    << abort(FatalError);
            }
        }
        else if (weightName != "none")
        {
            // Silently ignore "none", flag everything else as an error

            // TBD: treat missing "rho" like incompressible with rho = 1
            // and/or provided rhoRef value

            FatalErrorInFunction
                << "weightField " << weightName
                << " not found or an unsupported type" << nl
                << abort(FatalError);
        }
    }


    // Process the fields
    if (returnReduceOr(!vectorWeights.empty()))
    {
        if (scalarWeights.size())
        {
            vectorWeights *= scalarWeights;
        }

        writeAll(Sf, vectorWeights, points, faces);
    }
    else
    {
        writeAll(Sf, scalarWeights, points, faces);
    }


    if (operation_ != opNone)
    {
        file() << endl;
        Log << endl;
    }


    return true;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::updateMesh
(
    const mapPolyMesh& mpm
)
{
    needsUpdate_ = true;
}


void Foam::functionObjects::fieldValues::surfaceFieldValue::movePoints
(
    const polyMesh& mesh
)
{
    needsUpdate_ = true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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
#include "surfaceFields.H"
#include "polySurfaceFields.H"
#include "volFields.H"
#include "sampledSurface.H"
#include "surfaceWriter.H"
#include "interpolationCell.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class WeightType>
Foam::tmp<Foam::scalarField>
Foam::functionObjects::fieldValues::surfaceFieldValue::weightingFactor
(
    const Field<WeightType>& weightField,
    const bool useMag /* ignore */
)
{
    // The scalar form is specialized.
    // Other types: use mag() to generate a scalar field.
    return mag(weightField);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class WeightType>
inline bool Foam::functionObjects::fieldValues::surfaceFieldValue::canWeight
(
    const Field<WeightType>& fld
) const
{
    // Non-empty on some processor
    return returnReduceOr(!fld.empty());
}


template<class Type>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::validField
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;
    typedef DimensionedField<Type, polySurfaceGeoMesh> smt;

    return
    (
        foundObject<smt>(fieldName)
     || foundObject<vf>(fieldName)
     || (withSurfaceFields() && foundObject<sf>(fieldName))
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::getFieldValues
(
    const word& fieldName,
    const bool mandatory
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sf;
    typedef GeometricField<Type, fvPatchField, volMesh> vf;
    typedef DimensionedField<Type, polySurfaceGeoMesh> smt;

    if (foundObject<smt>(fieldName))
    {
        return lookupObject<smt>(fieldName);
    }
    else if (withSurfaceFields() && foundObject<sf>(fieldName))
    {
        return filterField(lookupObject<sf>(fieldName));
    }
    else if (foundObject<vf>(fieldName))
    {
        const vf& fld = lookupObject<vf>(fieldName);

        if (sampledPtr_)
        {
            // Could be runtime selectable
            // auto sampler = interpolation<Type>::New(sampleFaceScheme_, fld);

            // const interpolationCellPoint<Type> interp(fld);
            const interpolationCell<Type> interp(fld);

            return sampledPtr_->sample(interp);
        }
        else
        {
            return filterField(fld);
        }
    }

    if (mandatory)
    {
        FatalErrorInFunction
            << "Field " << fieldName << " not found in database" << nl
            << abort(FatalError);
    }

    return tmp<Field<Type>>::New();
}


template<class Type, class WeightType>
Type Foam::functionObjects::fieldValues::surfaceFieldValue::
processSameTypeValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const Field<WeightType>& weightField
) const
{
    Type result = Zero;
    switch (operation_)
    {
        case opNone:
        case typeScalar:
        case typeWeighted:
        case typeAbsolute:
        {
            break;
        }
        case opMin:
        {
            result = gMin(values);
            break;
        }
        case opMax:
        {
            result = gMax(values);
            break;
        }
        case opSumMag:
        {
            result = gSum(cmptMag(values));
            break;
        }
        case opSum:
        case opWeightedSum:
        case opAbsWeightedSum:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                tmp<scalarField> weight
                (
                    weightingFactor(weightField, Sf, is_magOp())
                );

                result = gSum(weight*values);
            }
            else
            {
                // Unweighted form
                result = gSum(values);
            }
            break;
        }
        case opSumDirection:
        case opSumDirectionBalance:
        {
            FatalErrorInFunction
                << "Operation " << operationTypeNames_[operation_]
                << " not available for values of type "
                << pTraits<Type>::typeName
                << exit(FatalError);

            break;
        }
        case opAverage:
        case opWeightedAverage:
        case opAbsWeightedAverage:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                const scalarField factor
                (
                    weightingFactor(weightField, Sf, is_magOp())
                );

                result = gSum(factor*values)/(gSum(factor) + ROOTVSMALL);
            }
            else
            {
                // Unweighted form
                const label n = returnReduce(values.size(), sumOp<label>());

                result = gSum(values)/(scalar(n) + ROOTVSMALL);
            }
            break;
        }
        case opAreaAverage:
        case opWeightedAreaAverage:
        case opAbsWeightedAreaAverage:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                const scalarField factor
                (
                    areaWeightingFactor(weightField, Sf, is_magOp())
                );

                result = gSum(factor*values)/gSum(factor + ROOTVSMALL);
            }
            else
            {
                // Unweighted form
                const scalarField factor(mag(Sf));

                result = gSum(factor*values)/gSum(factor + ROOTVSMALL);
            }
            break;
        }
        case opAreaIntegrate:
        case opWeightedAreaIntegrate:
        case opAbsWeightedAreaIntegrate:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                tmp<scalarField> factor
                (
                    areaWeightingFactor(weightField, Sf, is_magOp())
                );

                result = gSum(factor*values);
            }
            else
            {
                // Unweighted form
                tmp<scalarField> factor(mag(Sf));

                result = gSum(factor*values);
            }
            break;
        }
        case opCoV:
        {
            const scalarField magSf(mag(Sf));
            const scalar gSumMagSf = gSum(magSf);

            Type meanValue = gSum(values*magSf)/gSumMagSf;

            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                tmp<scalarField> vals(values.component(d));
                const scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res =
                    sqrt(gSum(magSf*sqr(vals - mean))/gSumMagSf)
                   /(mean + ROOTVSMALL);
            }

            break;
        }

        case opAreaNormalAverage:
        case opAreaNormalIntegrate:
        case opUniformity:
        {
            // Handled in specializations only
            break;
        }

        case opWeightedUniformity:
        case opAbsWeightedUniformity:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                // Change weighting from vector -> scalar and dispatch again
                return processValues<Type, scalar>
                (
                    values,
                    Sf,
                    weightingFactor(weightField, is_magOp())
                );
            }

            break;
        }
    }

    return result;
}


template<class Type, class WeightType>
Type Foam::functionObjects::fieldValues::surfaceFieldValue::processValues
(
    const Field<Type>& values,
    const vectorField& Sf,
    const Field<WeightType>& weightField
) const
{
    return processSameTypeValues(values, Sf, weightField);
}


template<class WeightType>
Foam::label Foam::functionObjects::fieldValues::surfaceFieldValue::writeAll
(
    const vectorField& Sf,
    const Field<WeightType>& weightField,
    const pointField& points,
    const faceList& faces
)
{
    label nProcessed = 0;

    // If using the surface writer, the points/faces parameters have already
    // been merged on the master and the writeValues routine will also gather
    // all data onto the master before calling the writer.
    // Thus only call the writer on master!

    // Begin writer time step
    if (Pstream::master() && surfaceWriterPtr_ && surfaceWriterPtr_->enabled())
    {
        auto& writer = *surfaceWriterPtr_;

        writer.open
        (
            points,
            faces,
            (
                outputDir()
              / regionTypeNames_[regionType_] + ("_" + regionName_)
            ),
            false  // serial - already merged
        );

        writer.beginTime(time_);
    }

    for (const word& fieldName : fields_)
    {
        if
        (
            writeValues<scalar>(fieldName, Sf, weightField, points, faces)
         || writeValues<vector>(fieldName, Sf, weightField, points, faces)
         || writeValues<sphericalTensor>
            (
                fieldName, Sf, weightField, points, faces
            )
         || writeValues<symmTensor>(fieldName, Sf, weightField, points, faces)
         || writeValues<tensor>(fieldName, Sf, weightField, points, faces)
        )
        {
            ++nProcessed;
        }
        else
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    // Finish writer time step
    if (Pstream::master() && surfaceWriterPtr_ && surfaceWriterPtr_->enabled())
    {
        auto& writer = *surfaceWriterPtr_;

        // Write geometry if no fields were written so that we still
        // can have something to look at.

        if (!writer.wroteData())
        {
            writer.write();
        }

        writer.endTime();
        writer.clear();
    }

    return nProcessed;
}


template<class Type, class WeightType>
bool Foam::functionObjects::fieldValues::surfaceFieldValue::writeValues
(
    const word& fieldName,
    const vectorField& Sf,
    const Field<WeightType>& weightField,
    const pointField& points,
    const faceList& faces
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(getFieldValues<Type>(fieldName, true));

        // Write raw values on surface if specified
        if (surfaceWriterPtr_ && surfaceWriterPtr_->enabled())
        {
            Field<Type> allValues(values);
            combineFields(allValues);

            if (Pstream::master())
            {
                fileName outputName =
                    surfaceWriterPtr_->write(fieldName, allValues);

                // Case-local file name with "<case>" to make relocatable
                dictionary propsDict;
                propsDict.add("file", time_.relativePath(outputName, true));
                this->setProperty(fieldName, propsDict);
            }
        }

        if (operation_ != opNone)
        {
            // Apply scale factor
            values *= scaleFactor_;

            Type result = processValues(values, Sf, weightField);

            switch (postOperation_)
            {
                case postOpSqrt:
                {
                    // sqrt: component-wise - does not change the type
                    for (direction d=0; d < pTraits<Type>::nComponents; ++d)
                    {
                        setComponent(result, d)
                            = sqrt(mag(component(result, d)));
                    }
                    break;
                }
                default:
                {
                    break;
                }
            }

            // Write state/results information
            word prefix, suffix;
            {
                if (postOperation_ != postOpNone)
                {
                    // Adjust result name to include post-operation
                    prefix += postOperationTypeNames_[postOperation_];
                    prefix += '(';
                    suffix += ')';
                }

                prefix += operationTypeNames_[operation_];
                prefix += '(';
                suffix += ')';
            }

            word resultName = prefix + regionName_ + ',' + fieldName + suffix;

            // Write state/results information

            Log << "    " << prefix << regionName_ << suffix
                << " of " << fieldName << " = ";


            // Operation or post-operation returns scalar?

            scalar sresult{0};

            bool alwaysScalar(operation_ & typeScalar);

            if (alwaysScalar)
            {
                sresult = component(result, 0);

                if (postOperation_ == postOpMag)
                {
                    sresult = mag(sresult);
                }
            }
            else if (postOperation_ == postOpMag)
            {
                sresult = mag(result);
                alwaysScalar = true;
            }


            if (alwaysScalar)
            {
                file()<< tab << sresult;

                Log << sresult << endl;

                this->setResult(resultName, sresult);
            }
            else
            {
                file()<< tab << result;

                Log << result << endl;

                this->setResult(resultName, result);
            }
        }
    }

    return ok;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    const labelList& own = field.mesh().faceOwner();
    const labelList& nei = field.mesh().faceNeighbour();

    auto tvalues = tmp<Field<Type>>::New(faceId_.size());
    auto& values = tvalues.ref();

    forAll(values, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];

        if (patchi >= 0)
        {
            // Boundary face - face id is the patch-local face id
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            // Internal face
            values[i] = 0.5*(field[own[facei]] + field[nei[facei]]);
        }
    }

    // No need to flip values - all boundary faces point outwards

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::surfaceFieldValue::filterField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    auto tvalues = tmp<Field<Type>>::New(faceId_.size());
    auto& values = tvalues.ref();

    forAll(values, i)
    {
        const label facei = faceId_[i];
        const label patchi = facePatchId_[i];

        if (patchi >= 0)
        {
            values[i] = field.boundaryField()[patchi][facei];
        }
        else
        {
            values[i] = field[facei];
        }
    }

    if (debug)
    {
        Pout<< "field " << field.name() << " oriented: "
            << field.is_oriented() << endl;
    }

    if (field.is_oriented())
    {
        forAll(values, i)
        {
            if (faceFlip_[i])
            {
                values[i] *= -1;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //

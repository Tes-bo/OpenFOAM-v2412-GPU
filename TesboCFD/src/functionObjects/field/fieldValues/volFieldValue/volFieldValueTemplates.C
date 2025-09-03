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

#include "volFieldValue.H"
#include "volFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::volFieldValue::validField
(
    const word& fieldName
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal IntVolFieldType;

    return
    (
        obr_.foundObject<VolFieldType>(fieldName)
     || obr_.foundObject<IntVolFieldType>(fieldName)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::fieldValues::volFieldValue::getFieldValues
(
    const word& fieldName,
    const bool mandatory
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal IntVolFieldType;

    if (obr_.foundObject<VolFieldType>(fieldName))
    {
        return filterField(obr_.lookupObject<VolFieldType>(fieldName));
    }
    else if (obr_.foundObject<IntVolFieldType>(fieldName))
    {
        return filterField(obr_.lookupObject<IntVolFieldType>(fieldName));
    }

    if (mandatory)
    {
        FatalErrorInFunction
            << "Field " << fieldName << " not found in database" << nl
            << abort(FatalError);
    }

    return tmp<Field<Type>>::New();
}


template<class Type>
Type Foam::functionObjects::fieldValues::volFieldValue::processValues
(
    const Field<Type>& values,
    const scalarField& V,
    const scalarField& weightField
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
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                result = gSum(weightField*values);
            }
            else
            {
                // Unweighted form
                result = gSum(values);
            }
            break;
        }
        case opAverage:
        case opWeightedAverage:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                result =
                    gSum(weightField*values)/(gSum(weightField) + ROOTVSMALL);
            }
            else
            {
                // Unweighted form
                const label n = returnReduce(values.size(), sumOp<label>());
                result = gSum(values)/(scalar(n) + ROOTVSMALL);
            }
            break;
        }
        case opVolAverage:
        case opWeightedVolAverage:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                result = gSum(weightField*V*values)
                    /(gSum(weightField*V) + ROOTVSMALL);
            }
            else
            {
                // Unweighted form
                result = gSum(V*values)/(gSum(V) + ROOTVSMALL);
            }
            break;
        }
        case opVolIntegrate:
        case opWeightedVolIntegrate:
        {
            if (is_weightedOp() && canWeight(weightField))
            {
                result = gSum(weightField*V*values);
            }
            else
            {
                // Unweighted form
                result = gSum(V*values);
            }
            break;
        }
        case opCoV:
        {
            const scalar sumV = gSum(V);

            Type meanValue = gSum(V*values)/sumV;

            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                tmp<scalarField> vals(values.component(d));
                const scalar mean = component(meanValue, d);
                scalar& res = setComponent(result, d);

                res = sqrt(gSum(V*sqr(vals - mean))/sumV)/(mean + ROOTVSMALL);
            }

            break;
        }
    }

    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::volFieldValue::writeValues
(
    const word& fieldName,
    const scalarField& V,
    const scalarField& weightField
)
{
    const bool ok = validField<Type>(fieldName);

    if (ok)
    {
        Field<Type> values(getFieldValues<Type>(fieldName));

        if (writeFields_)
        {
            word outName = fieldName + '_' + regionTypeNames_[regionType_];
            if (this->volRegion::regionName_ != polyMesh::defaultRegion)
            {
                outName = outName + '-' + this->volRegion::regionName_;
            }

            IOField<Type>
            (
                IOobject
                (
                    outName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                weightField.empty()
              ? scaleFactor_*values
              : scaleFactor_*weightField*values
            ).write();
        }

        if (operation_ != opNone)
        {
            // Apply scale factor
            values *= scaleFactor_;

            Type result = processValues(values, V, weightField);

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

            word regionPrefix;
            if (this->volRegion::regionName_ != polyMesh::defaultRegion)
            {
                regionPrefix = this->volRegion::regionName_ + ',';
            }

            word resultName = prefix + regionPrefix + fieldName + suffix;

            Log << "    " << prefix << this->volRegion::regionName_ << suffix
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
Foam::functionObjects::fieldValues::volFieldValue::filterField
(
    const Field<Type>& field
) const
{
    if (this->volRegion::useAllCells())
    {
        return field;
    }

    return tmp<Field<Type>>::New(field, cellIDs());
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(volFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, volFieldValue, runTime);
    addToRunTimeSelectionTable(functionObject, volFieldValue, dictionary);
}
}
}

const Foam::Enum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType
>
Foam::functionObjects::fieldValues::volFieldValue::operationTypeNames_
({
    // Normal operations
    { operationType::opNone, "none" },
    { operationType::opMin, "min" },
    { operationType::opMax, "max" },
    { operationType::opSum, "sum" },
    { operationType::opSumMag, "sumMag" },
    { operationType::opAverage, "average" },
    { operationType::opVolAverage, "volAverage" },
    { operationType::opVolIntegrate, "volIntegrate" },
    { operationType::opCoV, "CoV" },

    // Using weighting
    { operationType::opWeightedSum, "weightedSum" },
    { operationType::opWeightedAverage, "weightedAverage" },
    { operationType::opWeightedVolAverage, "weightedVolAverage" },
    { operationType::opWeightedVolIntegrate, "weightedVolIntegrate" },
});

const Foam::Enum
<
    Foam::functionObjects::fieldValues::volFieldValue::postOperationType
>
Foam::functionObjects::fieldValues::volFieldValue::postOperationTypeNames_
({
    { postOperationType::postOpNone, "none" },
    { postOperationType::postOpMag, "mag" },
    { postOperationType::postOpSqrt, "sqrt" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volFieldValue::usesVol()
const noexcept
{
    // Few operations require the cell volume
    switch (operation_)
    {
        case opVolAverage:
        case opVolIntegrate:
        case opWeightedVolAverage:
        case opWeightedVolIntegrate:
        case opCoV:
            return true;

        default:
            return false;
    }
}


void Foam::functionObjects::fieldValues::volFieldValue::writeFileHeader
(
    Ostream& os
) const
{
    volRegion::writeFileHeader(*this, os);

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

    // TBD: add in postOperation information?

    for (const word& fieldName : fields_)
    {
        os  << tab << operationTypeNames_[operation_]
            << "(" << fieldName << ")";
    }

    os  << endl;
}


Foam::label Foam::functionObjects::fieldValues::volFieldValue::writeAll
(
    const scalarField& V,
    const scalarField& weightField
)
{
    label nProcessed = 0;

    for (const word& fieldName : fields_)
    {
        if
        (
            writeValues<scalar>(fieldName, V, weightField)
         || writeValues<vector>(fieldName, V, weightField)
         || writeValues<sphericalTensor>(fieldName, V, weightField)
         || writeValues<symmTensor>(fieldName, V, weightField)
         || writeValues<tensor>(fieldName, V, weightField)
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

    return nProcessed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
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
    weightFieldNames_()
{
    read(dict);
    writeFileHeader(file());
}


Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
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
    weightFieldNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);

    weightFieldNames_.clear();

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

        Info<< nl << endl;
    }

    return true;
}


bool Foam::functionObjects::fieldValues::volFieldValue::write()
{
    volRegion::update();        // Ensure cached values are valid

    fieldValue::write();

    if (Pstream::master())
    {
        writeCurrentTime(file());
    }

    // Only some operations need the cell volume
    scalarField V;
    if (usesVol())
    {
        V = filterField(fieldValue::mesh_.V());
    }

    // Check availability and type of weight field
    // Only support a few weight types:
    // scalar: 0-N fields

    // Default is a zero-size scalar weight field (ie, weight = 1)
    scalarField scalarWeights;

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
    writeAll(V, scalarWeights);

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //

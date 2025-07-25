/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "fieldMinMax.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldMinMax::output
(
    const word& fieldName,
    const word& outputName,
    const label minCell,
    const label maxCell,
    const vector& minC,
    const vector& maxC,
    const label minProci,
    const label maxProci,
    const Type& minValue,
    const Type& maxValue
)
{
    OFstream& file = this->file();

    if (location_)
    {
        writeCurrentTime(file);

        writeTabbed(file, fieldName);

        file<< token::TAB << minValue
            << token::TAB << minC;

        if (Pstream::parRun())
        {
            file<< token::TAB << minProci;
        }

        file<< token::TAB << maxValue
            << token::TAB << maxC;

        if (Pstream::parRun())
        {
            file<< token::TAB << maxProci;
        }

        file<< endl;

        Log << "    min(" << outputName << ") = " << minValue
            << " in cell " << minCell
            << " at location " << minC;

        if (Pstream::parRun())
        {
            Log << " on processor " << minProci;
        }

        Log << nl << "    max(" << outputName << ") = " << maxValue
            << " in cell " << maxCell
            << " at location " << maxC;

        if (Pstream::parRun())
        {
            Log << " on processor " << maxProci;
        }
    }
    else
    {
        file<< token::TAB << minValue << token::TAB << maxValue;

        Log << "    min/max(" << outputName << ") = "
            << minValue << ' ' << maxValue;
    }

    Log << endl;

    // Write state/results information
    word nameStr('(' + outputName + ')');
    this->setResult("min" + nameStr, minValue);
    this->setResult("min" + nameStr + "_cell", minCell);
    this->setResult("min" + nameStr + "_position", minC);
    this->setResult("min" + nameStr + "_processor", minProci);
    this->setResult("max" + nameStr, maxValue);
    this->setResult("max" + nameStr + "_cell", maxCell);
    this->setResult("max" + nameStr + "_position", maxC);
    this->setResult("max" + nameStr + "_processor", maxProci);
}


template<class Type>
void Foam::functionObjects::fieldMinMax::calcMinMaxFieldType
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const word& outputFieldName
)
{
    const label proci = Pstream::myProcNo();

    // Find min/max internal field value info

    List<Type> minVs(Pstream::nProcs(), pTraits<Type>::max);
    labelList minCells(Pstream::nProcs(), Zero);
    List<vector> minCs(Pstream::nProcs(), Zero);

    List<Type> maxVs(Pstream::nProcs(), pTraits<Type>::min);
    labelList maxCells(Pstream::nProcs(), Zero);
    List<vector> maxCs(Pstream::nProcs(), Zero);

    labelPair minMaxIds = findMinMax(field);

    label minId = minMaxIds.first();
    if (minId != -1)
    {
        minVs[proci] = field[minId];
        minCells[proci] = minId;
        minCs[proci] = mesh_.C()[minId];
    }

    label maxId = minMaxIds.second();
    if (maxId != -1)
    {
        maxVs[proci] = field[maxId];
        maxCells[proci] = maxId;
        maxCs[proci] = mesh_.C()[maxId];
    }


    // Find min/max boundary field info
    const auto& fieldBoundary = field.boundaryField();
    const auto& CfBoundary = mesh_.C().boundaryField();

    forAll(fieldBoundary, patchi)
    {
        const Field<Type>& fp = fieldBoundary[patchi];
        if (fp.size())
        {
            const vectorField& Cfp = CfBoundary[patchi];

            const labelList& faceCells =
                fieldBoundary[patchi].patch().faceCells();

            minMaxIds = findMinMax(fp);

            minId = minMaxIds.first();
            if (minVs[proci] > fp[minId])
            {
                minVs[proci] = fp[minId];
                minCells[proci] = faceCells[minId];
                minCs[proci] = Cfp[minId];
            }

            maxId = minMaxIds.second();
            if (maxVs[proci] < fp[maxId])
            {
                maxVs[proci] = fp[maxId];
                maxCells[proci] = faceCells[maxId];
                maxCs[proci] = Cfp[maxId];
            }
        }
    }

    // Collect info from all processors and output
    Pstream::allGatherList(minVs);
    Pstream::allGatherList(minCells);
    Pstream::allGatherList(minCs);

    Pstream::allGatherList(maxVs);
    Pstream::allGatherList(maxCells);
    Pstream::allGatherList(maxCs);

    minId = findMin(minVs);
    const Type& minValue = minVs[minId];
    const label minCell = minCells[minId];
    const vector& minC = minCs[minId];

    maxId = findMax(maxVs);
    const Type& maxValue = maxVs[maxId];
    const label maxCell = maxCells[maxId];
    const vector& maxC = maxCs[maxId];

    output
    (
        field.name(),
        outputFieldName,
        minCell,
        maxCell,
        minC,
        maxC,
        minId,
        maxId,
        minValue,
        maxValue
    );
}


template<class Type>
void Foam::functionObjects::fieldMinMax::calcMinMaxFields
(
    const word& fieldName,
    const modeType& mode
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (obr_.foundObject<fieldType>(fieldName))
    {
        const fieldType& field = lookupObject<fieldType>(fieldName);

        switch (mode)
        {
            case mdMag:
            {
                calcMinMaxFieldType<scalar>
                (
                    mag(field),
                    word("mag(" + fieldName + ")")
                );
                break;
            }
            case mdCmpt:
            {
                calcMinMaxFieldType(field, fieldName);
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown min/max mode: " << modeTypeNames_[mode_]
                    << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "ccmReader.H"
#include "ccmInternal.H" // include last to avoid any strange interactions


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Get information about available fields.
// - assume that all fields are available for all solution intervals
// eg,
//   /FieldSets/FieldSet-1/Phase-1/Fields/Velocity
//   /FieldSets/FieldSet-1/Phase-1/Fields/Pressure
//   ...
void Foam::ccm::reader::determineFieldInfo
(
    const ccmID& fieldSetNode,
    fieldTable& table
)
{
    char fullName[kCCMIOMaxStringLength + 1];
    char shortName[kCCMIOProstarShortNameLength+1];

    // Loop through all phases
    ccmID phaseNode;
    int phaseI = 0;

    while
    (
        CCMIONextEntity
        (
            nullptr,
            fieldSetNode,
            kCCMIOFieldPhase,
            &phaseI,
            &phaseNode
        )
     == kCCMIONoErr
    )
    {
        CCMIODimensionality dims;

        // Examine each field in this fieldSet
        ccmID fieldNode;
        int fieldI = 0;

        // Get full/short names and dimension (scalar/vector/tensor).
        // Use short name as unique identifier.

        while
        (
            CCMIONextEntity
            (
                nullptr,
                phaseNode,
                kCCMIOField,
                &fieldI,
                &fieldNode
            )
         == kCCMIONoErr

         && CCMIOReadField
            (
                nullptr,
                fieldNode,
                fullName,
                shortName,
                &dims,
                nullptr
            )
         == kCCMIONoErr
        )
        {
            // The shortName *should* be a word, but some fields seem
            // to have blanks!
            {
                char *ptr = shortName;
                while (*ptr)
                {
                    if (!word::valid(*ptr))
                    {
                        *ptr = '_';
                    }
                    ptr++;
                }
            }

            word fieldName(shortName);

            if (dims == kCCMIOScalar)
            {
                // Add to table as required
                if (!table.found(fieldName))
                {
                    fieldEntry entry(fieldName, fullName);

                    entry.units
                    (
                        ccmReadOptstr("Units", fieldNode)
                    );

                    table.append(entry);
                }

                fieldEntry& entry = table.find(fieldName)();

                // Obtain sizes of data field
                // - only process cell/face data
                ccmID dataNode;
                int dataI = 0;

                CCMIODataLocation dataLocation;
                CCMIOIndex maxId;

                while
                (
                    CCMIONextEntity
                    (
                        nullptr,
                        fieldNode,
                        kCCMIOFieldData,
                        &dataI,
                        &dataNode
                    )
                 == kCCMIONoErr

                 && CCMIOEntitySize
                    (
                        nullptr,
                        dataNode,
                        nullptr,
                        &maxId
                    )
                 == kCCMIONoErr

                 && CCMIOReadFieldDatad
                    (
                        nullptr,
                        dataNode,
                        nullptr,
                        &dataLocation,
                        nullptr,
                        kCCMIOStart,
                        kCCMIOEnd
                    )
                 == kCCMIONoErr
                )
                {
                    if (dataLocation == kCCMIOCell)
                    {
                        entry.maxCellId(maxId);
                    }
                    else if (dataLocation == kCCMIOFace)
                    {
                        entry.maxFaceId(maxId);
                    }
                }
            }
        }
    }
}


// Get information about all available solutions.
// - it is sufficient to examine the first processor
//
// The "States":
//      Restart_1/Processor-1
//  -OR-
//      Transient_1/Processor-1
//      ...
//      Transient_N/Processor-1
//
// point to the "FieldSets":
//      FieldSet-1/RestartInfo
//
bool Foam::ccm::reader::detectSolution()
{
    // call once
    if (solutionStatus_ != UNKNOWN)
    {
        return (solutionStatus_ == OKAY || solutionStatus_ == READ);
    }

    // loop through all States
    ccmID stateNode;
    int stateI = 0;
    while
    (
        CCMIONextEntity
        (
            nullptr,
            (globalState_->root),
            kCCMIOState,
            &stateI,
            &stateNode
        )
     == kCCMIONoErr
    )
    {
        // Loop through all processors/solutions
        ccmID processorNode;
        ccmID solutionNode;
        int procI = 0;

        while
        (
            CCMIONextEntity
            (
                nullptr,
                stateNode,
                kCCMIOProcessor,
                &procI,
                &processorNode
            )
         == kCCMIONoErr

         && CCMIOReadProcessor
            (
                nullptr,
                processorNode,
                nullptr,        // Ignore verticesNode
                nullptr,        // Ignore topologyNode
                nullptr,        // Ignore initialField
                &solutionNode
            )
         == kCCMIONoErr
        )
        {
            // Restrict to solutions with RestartInfo on the first cpu
            // (there is normally only one set of restart data)
            ccmID restartNode;
            int restartI = 0;

            char  solutionName[kCCMIOMaxStringLength + 1];
            int   iteration = 0;
            float timeValue = 0;

            if
            (
                CCMIONextEntity
                (
                    nullptr,
                    solutionNode,
                    kCCMIORestart,
                    &restartI,
                    &restartNode
                )
             == kCCMIONoErr

             && CCMIOEntityName
                (
                    nullptr,
                    stateNode,
                    solutionName
                )
             == kCCMIONoErr

             && CCMIOReadRestartInfo
                (
                    nullptr,
                    restartNode,
                    nullptr,       // Ignore solverName
                    &iteration,
                    &timeValue,
                    nullptr,       // Ignore timeUnits
                    nullptr        // Ignore startAngle
                )
             == kCCMIONoErr
            )
            {
                solutionTable_.append
                (
                    solutionEntry(solutionName, iteration, timeValue)
                );
            }

            // Determine field information
            determineFieldInfo(solutionNode, fieldTable_);

            // check Lagrangian data
            ccmID lagrangianNode;
            ccmID lagrangianSolutions;
            int lagrangianI = 0;

            if
            (
                CCMIONextEntity
                (
                    nullptr,
                    processorNode,
                    kCCMIOLagrangianData,
                    &lagrangianI,
                    &lagrangianNode
                )
             == kCCMIONoErr

             && CCMIOReadLagrangianData
                (
                    nullptr,
                    lagrangianNode,
                    nullptr,
                    &lagrangianSolutions
                )
             == kCCMIONoErr
            )
            {
                determineFieldInfo(lagrangianSolutions, lagrangianTable_);
            }
        }
    }

    if (solutionTable_.size() && fieldTable_.size())
    {
        solutionStatus_ = OKAY;
    }
    else
    {
        solutionStatus_ = BAD;
    }

    return (solutionStatus_ == OKAY || solutionStatus_ == READ);
}

#define SOLID_STRESS_HACK

//
// Read solution and field combination
//
Foam::tmp<Foam::scalarField>
Foam::ccm::reader::readField
(
    const word& solutionName,
    const word& fieldName,
    const bool wallData
)
{
    // get State by name
    ccmID stateNode;
    ccmID solutionNode;

    if
    (
        CCMIOGetState
        (
            nullptr,
            (globalState_->root),
            solutionName.c_str(),
            nullptr,
            &stateNode
        )
     != kCCMIONoErr

     || !fieldTable_.found(fieldName)
    )
    {
        return tmp<scalarField>::New();
    }

    CCMIODataLocation requestedLocation = kCCMIOCell;

    fieldEntry& entry = fieldTable_.find(fieldName)();

    label maxId = entry.maxCellId();

    if (wallData)
    {
        maxId = entry.maxFaceId();
        requestedLocation = kCCMIOFace;
    }

    // we can skip empty fields immediately
    if (!maxId)
    {
        return tmp<scalarField>::New();
    }

    char shortName[kCCMIOProstarShortNameLength+1];

    List<label>  mapData;
    List<scalar> rawData;

    tmp<scalarField> tscalarData
    (
        new scalarField(maxId + 1, option().undefScalar())
    );
    scalarField& scalarData = tscalarData.ref();


    CCMIODimensionality dims;

    // Loop across all processors
    ccmID processorNode;
    int procI = 0;

    while
    (
        CCMIONextEntity
        (
            nullptr,
            stateNode,
            kCCMIOProcessor,
            &procI,
            &processorNode
        )
     == kCCMIONoErr

    && CCMIOReadProcessor
        (
            nullptr,
            processorNode,
            nullptr,        // Ignore verticesNode
            nullptr,        // Ignore topologyNode
            nullptr,        // Ignore initialField
            &solutionNode
        )
     == kCCMIONoErr
    )
    {
        // loop through all phases
        int phaseI = 0;
        ccmID phaseNode;

        while
        (
            CCMIONextEntity
            (
                nullptr,
                solutionNode,
                kCCMIOFieldPhase,
                &phaseI,
                &phaseNode
            )
         == kCCMIONoErr
        )
        {
            // Get Field that matches the name
            ccmID fieldNode;
            int fieldI = 0;

            while
            (
                CCMIONextEntity
                (
                    nullptr,
                    phaseNode,
                    kCCMIOField,
                    &fieldI,
                    &fieldNode
                )
                == kCCMIONoErr
            )
            {
                // Get full/short names and dimension (scalar/vector/tensor)
                // use short name as unique identifier
                CCMIOReadField
                (
                    &(globalState_->error),
                    fieldNode,
                    nullptr,
                    shortName,
                    &dims,
                    nullptr
                );
                assertNoError
                (
                    "reading post data field: "
                  + string(shortName)
                );

                if (fieldName == shortName && dims == kCCMIOScalar)
                {
                    // Obtain sizes of data field
                    ccmID dataNode;
                    ccmID mapId;
                    int dataI = 0;
                    while
                    (
                        CCMIONextEntity
                        (
                            nullptr,
                            fieldNode,
                            kCCMIOFieldData,
                            &dataI,
                            &dataNode
                        )
                     == kCCMIONoErr
                    )
                    {
                        CCMIODataLocation dataLocation;
                        CCMIOSize  n;

                        // Only process cell/face data
                        if
                        (
                            CCMIOEntitySize
                            (
                                nullptr,
                                dataNode,
                                &n,
                                nullptr
                            )
                         == kCCMIONoErr

                         && CCMIOReadFieldDatad
                            (
                                nullptr,
                                dataNode,
                                &mapId,
                                &dataLocation,
                                nullptr,
                                kCCMIOStart,
                                kCCMIOEnd
                            )
                         == kCCMIONoErr

                         && dataLocation == requestedLocation
                        )
                        {
#ifdef SOLID_STRESS_HACK
                            bool okayCombination = true;

                            CCMIOSize len;
                            if
                            (
                                CCMIOEntityDescription
                                (
                                    nullptr,
                                    dataNode,
                                    &len,
                                    nullptr
                                )
                             == kCCMIONoErr
                            )
                            {
                                char* dataLabel = new char[len + 1];

                                if
                                (
                                    CCMIOEntityDescription
                                    (
                                        nullptr,
                                        dataNode,
                                        &len,
                                        dataLabel
                                    )
                                 == kCCMIONoErr
                                )
                                {
                                    if
                                    (
                                        (
                                            strstr(fieldName.c_str(), "SIG")
                                         || strstr(fieldName.c_str(), "EPS")
                                        )
                                     && strstr(dataLabel,  "So") == nullptr
                                    )
                                    {
                                        okayCombination = false;
                                        // skip non-solid
                                    }
                                }

                                delete[] dataLabel;
                            }

                            if (!okayCombination)
                            {
                                continue;
                            }
#endif

                            mapData.setSize(n);
                            rawData.setSize(n);

                            readMap
                            (
                                mapId,
                                mapData
                            );

                            CCMIOReadFieldDatad
                            (
                                &(globalState_->error),
                                dataNode,
                                nullptr,
                                nullptr,
                                rawData.data(),
                                kCCMIOStart,
                                kCCMIOEnd
                            );
                            assertNoError
                            (
                                "reading post data field: "
                              + string(shortName)
                            );

                            // transcribe to output list
                            forAll(mapData, i)
                            {
                                const label cellId = mapData[i];
                                scalarData[cellId] = rawData[i];
                            }
                        }
                    }
                }
            }
        }
    }

    // Overwrite possible junk in cell 0 (often used as /dev/null)
    scalarData[0] = option().undefScalar();

    return tscalarData;
}


// ************************************************************************* //

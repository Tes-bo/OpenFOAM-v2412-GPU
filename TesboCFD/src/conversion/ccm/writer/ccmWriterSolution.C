/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "ccmWriter.H"
#include "dictionary.H"
#include "IFstream.H"
#include "StringStream.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ccmInternal.H"  // include last to avoid any strange interactions


// * * * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * * //

// names for known field types
Foam::dictionary Foam::ccm::writer::defaultNameMapping
(
    IStringStream
    (
        // volScalarField
        "p   { name P; Label \"Pressure\"; units \"Pa\"; }"
     // "?? { name PTER; Label \"Thermo. Pressure\"; units \"Pa\"; }"
     // "yPlus { name YPLU; Label \"YPLUS\"; }"
     // "?? { name DIST; Label \"Normal Distance\"; units \"m\"; }"
     // "?? { name H; Label \"Enthalpy\"; units \"J/kg\"; }"
        "rho { name DENS; Label \"Density\"; units \"kg/m3\"; }"
        "T   { name T; Label \"Temperature\"; units \"K\"; }"
        "k   { name TE; Label \"Turb Kinetic Energy\"; units \"m2/s2\"; }"
        "epsilon { name ED; Label \"Dissipation Rate\"; units \"m2/s3\"; }"
        "mu  { name LAMV; Label \"Molecular Viscosity\"; units \"Pa s\"; }"
        "mut { name VIS;  Label \"Turbulent Viscosity\"; units \"Pa s\"; }"
        // volVectorField
        "U   { name ALL; Label \"Velocity\"; units \"m/s\"; }"
        // U-components:
        "_0U { name SU; Label \"Velocity Component U\"; units \"m/s\"; }"
        "_1U { name SV; Label \"Velocity Component V\"; units \"m/s\"; }"
        "_2U { name SW; Label \"Velocity Component W\"; units \"m/s\"; }"
        // surfaceScalarField
        "phi { name MassFlux; Label \"Mass Flux\"; units \"kg/s\"; }"
    )()
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::ccm::writer::newFieldNode
(
    const ccmID& phaseNode,
    const word& fieldName,
    const dictionary& nameMapping,
    const ccmDimension& dims,
    ccmID& fieldNode
) const
{
    const dictionary* subDictPtr = nameMapping.findDict(fieldName);

    if (!subDictPtr)
    {
        return false;
    }

    const dictionary& dict = *subDictPtr;

    word shortName;
    if (!dict.readIfPresent("name", shortName))
    {
        return false;
    }

    string ccmLabel(fieldName);
    dict.readIfPresent("Label", ccmLabel);

    CCMIONewField
    (
        &(globalState_->error),
        phaseNode,
        ccmLabel.c_str(),
        shortName.c_str(),
        dims(),
        &fieldNode
    );


    string units;
    if (dims() == kCCMIOScalar && dict.readIfPresent("units", units))
    {
        CCMIOWriteOptstr
        (
            &(globalState_->error),
            fieldNode,
            "Units",
            units.c_str()
        );
    }

    return true;

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ccm::writer::writeSolution
(
    const IOobjectList& objects,
    const fileName& remappingDictName
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if (!isA<fvMesh>(mesh_))
    {
        WarningInFunction
            << "cannot write solutions with a polyMesh instead of a fvMesh"
            << endl;
        return;
    }

    dictionary remapDict;

    if (remappingDictName.empty())
    {
        remapDict = IOdictionary
        (
            IOobject
            (
                "remapping",
                mesh_.time().constant(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
    }
    else
    {
        // Specified (absolute/relative) name: treat like MUST_READ
        remapDict = dictionary(IFstream(remappingDictName)());
    }


    dictionary nameMapping(defaultNameMapping);

    // Merge without overwrite
    const dictionary* subDictPtr = remapDict.findDict("fields");

    if (subDictPtr)
    {
        nameMapping |= *subDictPtr;
    }


    ccmID stateNode, processorNode, nodeId;
    int procI = 0;

    // Create first ("default") state node
    CCMIONewState
    (
        &(globalState_->error),
        (globalState_->root),
        "default",
        nullptr,
        nullptr,
        &stateNode
    );
    assertNoError("could not create default state");

    // Use first processor or create a new one
    procI = 0;
    if
    (
        CCMIONextEntity
        (
            nullptr,
            stateNode,
            kCCMIOProcessor,
            &procI,
            &processorNode
        )
     != kCCMIONoErr
    )
    {
        CCMIONewEntity
        (
            &(globalState_->error),
            stateNode,
            kCCMIOProcessor,
            nullptr,
            &processorNode
        );
        assertNoError("could not create processor node");
    }

    // Remove old data (if it exists)
    CCMIOClearProcessor
    (
        nullptr, stateNode, processorNode,
        true,   // Clear vertices
        true,   // Clear topology
        true,   // Clear initial field
        true,   // Clear solution
        true    // Clear lagrangian
    );

    //
    // Create vertices and topology nodes
    // and references to vertices and topology files
    //
    {
        ccmID verticesNode, topoNode;
        CCMIONewEntity
        (
            &(globalState_->error),
            (globalState_->root),
            kCCMIOVertices,
            nullptr,
            &verticesNode
        );

        CCMIONewEntity
        (
            &(globalState_->error),
            (globalState_->root),
            kCCMIOTopology,
            nullptr,
            &topoNode
        );

        string topoFileName = defaultMeshName + ".ccmg";

        CCMIOWriteProcessor
        (
            &(globalState_->error),
            processorNode,
            topoFileName.c_str(), &verticesNode,// verticesFile, verticesNode
            topoFileName.c_str(), &topoNode,    // topologyFile, topologyNode
            nullptr, nullptr,   // initialField unchanged
            nullptr, nullptr    // no solutionFile, solution unchanged
        );
    }
    assertNoError("Error after writing geometry processor");

    CCMIONewState
    (
        &(globalState_->error),
        (globalState_->root),
        "Restart_1",
        nullptr,
        nullptr,
        &stateNode
    );
    assertNoError("could not create Restart_1 state");

    // Use first processor or create a new one
    procI = 0;
    if
    (
        CCMIONextEntity
        (
            nullptr,
            stateNode,
            kCCMIOProcessor,
            &procI,
            &processorNode
        )
     != kCCMIONoErr
    )
    {
        CCMIONewEntity
        (
            &(globalState_->error),
            stateNode,
            kCCMIOProcessor,
            nullptr,
            &processorNode
        );
        assertNoError("could not create 'Processor' node");
    }

    // Remove old data (if it exists)
    CCMIOClearProcessor
    (
        nullptr, stateNode, processorNode,
        true,   // Clear vertices
        true,   // Clear topology
        true,   // Clear initial field
        true,   // Clear solution
        true    // Clear lagrangian
    );

    // Write out some simple solution data
    ccmID phaseNode, fieldSetNode;

    // Use first FieldSet
    int fieldSetI = 0;
    if
    (
        CCMIONextEntity
        (
            nullptr,
            (globalState_->root),
            kCCMIOFieldSet,
            &fieldSetI,
            &fieldSetNode
        )
     != kCCMIONoErr
    )
    {
        CCMIONewEntity
        (
            &(globalState_->error),
            (globalState_->root),
            kCCMIOFieldSet,
            nullptr,
            &fieldSetNode
        );
        assertNoError("could not create FieldSet node");
    }

    // RestartInfo
    {
        ccmID restartInfoNode, restartDataNode;
        CCMIONewEntity
        (
            &(globalState_->error),
            fieldSetNode,
            kCCMIORestart,
            nullptr,
            &restartInfoNode
        );
        assertNoError("could not create restartInfoNode node");

        // Get time information
        const Time& runTime = mesh_.time();
        label timeIndex = 0;
        if
        (
            runTime.timeName() != runTime.constant()
         && runTime.timeName() != "0"
        )
        {
            IOobject io
            (
                "time",
                runTime.timeName(),
                "uniform",
                runTime,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            );

            if (io.typeHeaderOk<IOdictionary>(true))
            {
                IOdictionary(io).readEntry("index", timeIndex);
            }
        }

        CCMIOWriteRestartInfo
        (
            &(globalState_->error),
            restartInfoNode,
            "ccm::writer",      // solverName
            timeIndex,          // iteration
            0.0,                // time
            nullptr,            // timeUnits: default  (s)
            0.0                 // startAngle: non-rotating mesh
        );

        // RestartData
        CCMIONewEntity
        (
            &(globalState_->error),
            restartInfoNode,
            kCCMIORestartData,
            nullptr,
            &restartDataNode
        );
        assertNoError("could not create restartDataNode node");

        // Minimal information required by PROSTAR

        CCMIOWriteOptf
        (
            nullptr,
            restartDataNode,
            "SCALE",
            0.001        // [m] -> [mm]: PROSTAR is still a bit buggy here
        );


        // Add Phase data
        CCMIONewIndexedEntity
        (
            &(globalState_->error),
            restartDataNode,
            kCCMIOFieldPhase,
            1,
            nullptr,
            &nodeId
        );

        // HACK:  for calculating Mach number
        // FIXME: use thermodynamicProperties /or/ thermophysicalProperties
        CCMIOWriteOptf
        (
            nullptr,
            nodeId,
            "Material Specific Heat",
            1007
        );

        CCMIOWriteOptf
        (
            nullptr,
            nodeId,
            "Material Molecular Weight",
            28.96
        );
    }

    CCMIONewIndexedEntity
    (
        &(globalState_->error),
        fieldSetNode,
        kCCMIOFieldPhase,
        1,
        "Phase_0001",
        &phaseNode
    );

    forAllConstIters(objects, iter)
    {
        const word fieldName = (*iter()).name();
        const word fieldType = (*iter()).headerClassName();

        if (!nameMapping.found(fieldName))
        {
            continue;
        }

        if (fieldType == volScalarField::typeName)
        {
            Info<< " " << fieldName << flush;

            volScalarField field
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                refCast<const fvMesh>(mesh_)
            );

            ccmID fieldNode;

            // Info<< " call newFieldNode " << fieldName << flush;

            if
            (
                !newFieldNode
                (
                    phaseNode,
                    fieldName,
                    nameMapping,
                    kCCMIOScalar,
                    fieldNode
                )
            )
            {
                continue;
            }

            // Write cellData first
            CCMIONewEntity
            (
                &(globalState_->error),
                fieldNode,
                kCCMIOFieldData,
                nullptr,
                &nodeId
            );

            CCMIOWriteFieldDatad
            (
                &(globalState_->error),
                nodeId,
                maps_->cells,
                kCCMIOCell,
                const_cast<scalar*>
                (
                    field.primitiveField().begin()
                ),
                kCCMIOStart,
                kCCMIOEnd
            );
            assertNoError("writing internalField " + fieldName);

            // Write boundaryData
            forAll(patches, patchI)
            {
                CCMIONewEntity
                (
                    &(globalState_->error),
                    fieldNode,
                    kCCMIOFieldData,
                    nullptr,
                    &nodeId
                );

                CCMIOWriteFieldDatad
                (
                    &(globalState_->error),
                    nodeId,
                    maps_->boundary[patchI],
                    kCCMIOFace,
                    const_cast<scalar*>
                    (
                        field.boundaryField()[patchI].begin()
                    ),
                    kCCMIOStart,
                    kCCMIOEnd
                );
            }

            assertNoError("writing boundaryField " + fieldName);
        }
        else if (fieldType == volVectorField::typeName && fieldName == "U")
        {
            Info<< " " << fieldName << flush;

            volVectorField vfield
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                refCast<const fvMesh>(mesh_)
            );


            ccmID vectorNode;
            newFieldNode
            (
                phaseNode,
                fieldName,
                nameMapping,
                kCCMIOVector,
                vectorNode
            );

            for (direction cmpt=0; cmpt < pTraits<vector>::nComponents; ++cmpt)
            {
                word componentName("_" + Foam::name(cmpt) + fieldName);
                volScalarField field(vfield.component(cmpt));

                CCMIOComponent ccmComponent = kCCMIOVectorX;
                switch(cmpt) {
                    case 0:
                        ccmComponent = kCCMIOVectorX;
                        break;
                    case 1:
                        ccmComponent = kCCMIOVectorY;
                        break;
                    case 2:
                        ccmComponent = kCCMIOVectorZ;
                        break;
                }

                ccmID componentNode;
                if
                (
                    !newFieldNode
                    (
                        phaseNode,
                        componentName,
                        nameMapping,
                        kCCMIOScalar,
                        componentNode
                    )
                )
                {
                    continue;
                }

                // Re-register with vector field
                CCMIOWriteMultiDimensionalFieldData
                (
                    &(globalState_->error),
                    vectorNode,
                    ccmComponent,
                    componentNode
                );

                // Write cellData first
                CCMIONewEntity
                (
                    &(globalState_->error),
                    componentNode,
                    kCCMIOFieldData,
                    nullptr,
                    &nodeId
                );

                CCMIOWriteFieldDatad
                (
                    &(globalState_->error),
                    nodeId,
                    maps_->cells,
                    kCCMIOCell,
                    const_cast<scalar*>
                    (
                        field.primitiveField().begin()
                    ),
                    kCCMIOStart, kCCMIOEnd
                );
                assertNoError
                (
                    "writing internalField " + fieldName + " " + componentName
                );


                // Write boundaryData
                forAll(patches, patchI)
                {
                    CCMIONewEntity
                    (
                        &(globalState_->error),
                        componentNode,
                        kCCMIOFieldData,
                        nullptr,
                        &nodeId
                    );

                    CCMIOWriteFieldDatad
                    (
                        &(globalState_->error),
                        nodeId,
                        maps_->boundary[patchI],
                        kCCMIOFace,
                        const_cast<scalar*>
                        (
                            field.boundaryField()[patchI].begin()
                        ),
                        kCCMIOStart, kCCMIOEnd
                    );
                }

                assertNoError
                (
                    "writing boundaryField " + fieldName + " " + componentName
                );
            }
        }
        else if (fieldType == surfaceScalarField::typeName)
        {
#if 0
            // Still have problems reading surface fields in PROSTAR
            Info<< " " << fieldName << flush;

            surfaceScalarField field
            (
                IOobject
                (
                    fieldName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                refCast<const fvMesh>(mesh_)
            );

            ccmID fieldNode;

            // Info<< " call newFieldNode " << fieldName << flush;

            if
            (
                !newFieldNode
                (
                    phaseNode,
                    fieldName,
                    nameMapping,
                    kCCMIOScalar,
                    fieldNode
                )
            )
            {
                continue;
            }

            // Write cell faceData first
            CCMIONewEntity
            (
                &(globalState_->error),
                fieldNode,
                kCCMIOFieldData,
                nullptr,
                &nodeId
            );

            CCMIOWriteFieldDatad
            (
                &(globalState_->error),
                nodeId,
                maps_->internalFaces,
                kCCMIOFace,
                const_cast<scalar*>
                (
                    field.primitiveField().cdata()
                ),
                kCCMIOStart,
                kCCMIOEnd
            );
            assertNoError("writing internalField " + fieldName);

            // Write boundaryData
            forAll(patches, patchI)
            {
                CCMIONewEntity
                (
                    &(globalState_->error),
                    fieldNode,
                    kCCMIOFieldData,
                    nullptr,
                    &nodeId
                );

                CCMIOWriteFieldDatad
                (
                    &(globalState_->error),
                    nodeId,
                    maps_->boundary[patchI],
                    kCCMIOFace,
                    const_cast<scalar*>
                    (
                        field.boundaryField()[patchI].cdata()
                    ),
                    kCCMIOStart,
                    kCCMIOEnd
                );
            }

            assertNoError("writing boundaryField " + fieldName);
#endif
        }
    }
    Info<< endl;

    assertNoError("Error before writing processor");

    CCMIOWriteProcessor
    (
        &(globalState_->error),
        processorNode,
        nullptr, nullptr,     // vertices
        nullptr, nullptr,     // topology
        nullptr, nullptr,     // initial field
        nullptr, &fieldSetNode
    );
    assertNoError("Error after writing processor");
}


// ************************************************************************* //

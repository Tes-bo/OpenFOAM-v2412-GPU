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
#include "IFstream.H"
#include "IOdictionary.H"

#include "ccmInternal.H"  // include last to avoid any strange interactions

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The cellTable integer options
const char* Foam::ccm::reader::cellTableOpti[] =
{
    "MaterialId", "PorosityId", "SpinId", "GroupId", "ColorIdx", nullptr
};

// The cellTable string options - "Label" handled separately
const char* Foam::ccm::reader::cellTableOptstr[] =
{
    "MaterialType", nullptr
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Simulate GetEntityIndex
// CCMIOGetEntityIndex(nullptr, nodeId, &internalId);
int Foam::ccm::reader::ccmGetEntityIndex(ccmNODE node)
{
    int intval = 0;
    char name[kCCMIOMaxStringLength + 1];

    if
    (
        CCMIOGetName
        (
            nullptr,
            node,
            name
        )
        == kCCMIONoErr
    )
    {
        char const *pos = strrchr(name, '-');

        if (!pos)
        {
            intval = 0;
        }
        else
        {
            // Negative ID number
            if (pos != name && (*(pos - 1) == '-'))
                --pos;

            intval = atoi(++pos);
        }
    }

    return intval;
}


std::string Foam::ccm::reader::ccmReadNodestr
(
    const char* opt,
    ccmNODE node
)
{
    char *strval = 0;
    std::string str;

    if
    (
        CCMIOReadNodestr
        (
            nullptr,
            node,
            opt,
            &strval
        )
     == kCCMIONoErr
     && strval
    )
    {
        str = strval;
    }

    // Allocated by CCMIOReadNodestr
    if (strval)
    {
        free(strval);
    }

    return str;
}


std::string Foam::ccm::reader::ccmReadOptstr
(
    const char* opt,
    ccmID node
)
{
    std::string str;
    int len = 0;

    if
    (
        CCMIOReadOptstr
        (
            nullptr,
            node,
            opt,
            &len,
            nullptr
        )
     == kCCMIONoErr
     && len
    )
    {
        char* strval = new char[len + 1];
        if
        (
            CCMIOReadOptstr
            (
                nullptr,
                node,
                opt,
                &len,
                strval
            )
         == kCCMIONoErr
        )
        {
            strval[len] = '\0';
            str = strval;
        }
        delete[] strval;
    }

    return str;
}


// Read map data and check error
void Foam::ccm::reader::readMap
(
    const ccmID& mapId,
    labelList& data
)
{
    if (globalState_->hasError())
    {
        return;
    }

    CCMIOReadMap
    (
        &(globalState_->error),
        mapId,
        data.data(),
        kCCMIOStart,
        kCCMIOEnd
    );
    assertNoError("error reading map");
}


// readProblemDescription:
// - get cellTable and boundaryRegion information
// - some regions defined here may be unused by our mesh
void Foam::ccm::reader::readProblemDescription
(
    const ccmID& probNode
)
{
    readInterfaceDefinitions();
    readProblemDescription_boundaryRegion(probNode);
    readProblemDescription_cellTable(probNode);

#ifdef DEBUG_CCMIOREAD
    Info<< "InterfaceDefinitions: " << interfaceDefinitions_ << nl
        << "cellTable" << cellTable_ << nl
        << "boundaryRegion" << boundaryRegion_ << endl;
#endif
}


// readInterfaceDefinitions:
// - get /InterfaceDefinitions.
//   used by STARCCM to define in-place interfaces, etc
// - only handle in-place one here
void Foam::ccm::reader::readInterfaceDefinitions()
{
    interfaceDefinitions_.clear();

    // /InterfaceDefinitions
    // ---------------------

    CCMIONode rootNode;
    CCMIONode interfaceDefNode;

    // CCMIOID -> CCMIONODE
    if
    (
        CCMIOGetEntityNode
        (
            nullptr,
            (globalState_->root),
            &rootNode
        )
     == kCCMIONoErr

        // Get "/InterfaceDefinitions"
     &&
        CCMIOGetNode
        (
            nullptr,
            rootNode,
            "InterfaceDefinitions",
            &interfaceDefNode
        )
     == kCCMIONoErr
    )
    {
        CCMIONode interfaceNode;

        // Simulate CCMIONextEntity
        for
        (
            int index = 0;
            CCMIOGetNextChildWithLabel
            (
                nullptr,
                interfaceDefNode,
                "Interface",
                &index,
                &interfaceNode
            ) == kCCMIONoErr;
            /* nop */
        )
        {
            interfaceEntry ifentry(ccmGetEntityIndex(interfaceNode));

            if
            (
                CCMIOReadNodei
                (
                    nullptr,
                    interfaceNode,
                    "Boundary0",
                    &(ifentry.bnd0)
                )
             == kCCMIONoErr

             && CCMIOReadNodei
                (
                    nullptr,
                    interfaceNode,
                    "Boundary1",
                    &(ifentry.bnd1)
                )
             == kCCMIONoErr

                // only handle 'IN_PLACE' interfaces
             && interfaceEntry::isInPlace
                (
                    ccmReadNodestr("Configuration", interfaceNode)
                )
            )
            {
                interfaceDefinitions_.add(ifentry);
            }
        }
    }
}


// readProblemDescription - get boundaryRegion information
//
// CCM Node: /ProblemDescriptions/ProblemDescription-N/BoundaryRegion-N
//
// Requires InterfaceDefinitions first
void Foam::ccm::reader::readProblemDescription_boundaryRegion
(
    const ccmID& probNode
)
{
    if (option().useNumberedNames())
    {
        Info<< "using numbered patch/zone names" << endl;
    }

    boundaryRegion_.clear();

    ccmID node;

    // boundaryRegion information
    // --------------------------
    for
    (
        int nodeI = 0;
        CCMIONextEntity
        (
            nullptr,
            probNode,
            kCCMIOBoundaryRegion,
            &nodeI,
            &node
        ) == kCCMIONoErr;
        /* nop */
    )
    {
        // Read boundaryRegionId
        int Id = 0;
        CCMIOGetEntityIndex
        (
            &(globalState_->error),
            node,
            &Id
        );
        assertNoError("error reading boundaryRegion index");

        dictionary dict;

        // Read boundary type
        // (inlet|outlet|symmetry|wall|
        //   cyclic|stagnation|pressure|baffle|freestream)
        {
            const char* opt = "BoundaryType";
            std::string str = ccmReadOptstr(opt, node);

            if (str.empty())
            {
                dict.add(opt, "empty");
            }
            else if (str == "internal")
            {
                // Old PROSTAR bug: "monitoring" mislabeled as "internal"
                dict.add(opt, "monitoring");
            }
            else
            {
                dict.add(opt, word::validate(str, true));
            }
        }


        // Read boundary name:
        // - from 'Label' field (STARCCM+)
        // - from 'BoundaryName' field (PROSTAR and/or STARCCM+)
        // - fallback: generate one from boundary type and boundaryRegionId
        {
            const char* opt = "Label";
            std::string str;

            if (option().useNumberedNames())
            {
                str = "patch_" + ::Foam::name(Id);
            }
            else if
            (
                option().renameInterfaces()
             && interfaceDefinitions_.isInterface(Id)
            )
            {
#ifdef DEBUG_CCMIOREAD
                Info<< "boundary is on an interface: remap name for  "
                    << Id << endl;
#endif
                // Substitute immediately with interface name
                str = interfaceDefinitions_.interfaceName(Id);
            }
            else if
            (
                (str = ccmReadOptstr(opt, node)).empty()
             && (str = ccmReadOptstr("BoundaryName", node)).empty()
            )
            {
                // Fallback
                str = dict.get<word>("BoundaryType") + "_" + ::Foam::name(Id);
            }

            if (!str.empty())
            {
                dict.add(opt, word::validate(str, true));
            }
        }

        boundaryRegion_.insert(Id, dict);
    }
}


// readProblemDescription - get cellTable information
//
// CCM Node: /ProblemDescriptions/ProblemDescription-N/CellType-N
//
void Foam::ccm::reader::readProblemDescription_cellTable
(
    const ccmID& probNode
)
{
    cellTable_.clear();

    ccmID node;

    // cellTable information
    // ---------------------
    for
    (
        int nodeI = 0;
        CCMIONextEntity
        (
            nullptr,
            probNode,
            kCCMIOCellType,
            &nodeI,
            &node
        ) == kCCMIONoErr;
        /* nop */
    )
    {
        // Read cellTableId (CellType)
        int Id = 0;
        CCMIOGetEntityIndex
        (
            &(globalState_->error),
            node,
            &Id
        );
        assertNoError("error reading cellTable index");

        dictionary dict;

        // Special treatment for "Label"
        {
            const char* opt = "Label";
            std::string str;

            if (!option().useNumberedNames())
            {
                str = ccmReadOptstr(opt, node);
            }

            if (str.empty())
            {
                str = "zone_" + ::Foam::name(Id);
            }

            dict.add(opt, word::validate(str, true));
        }


        // Other string options
        for (int i=0; cellTableOptstr[i]; ++i)
        {
            const char* opt = cellTableOptstr[i];
            std::string str = ccmReadOptstr(opt, node);

            if (!str.empty())
            {
                dict.add(opt, word::validate(str, true));
            }
        }

        // Add non-zero integer options
        for (int i=0; cellTableOpti[i]; ++i)
        {
            const char* opt = cellTableOpti[i];
            int intval;

            if
            (
                CCMIOReadOpti
                (
                    nullptr,
                    node,
                    opt,
                    &intval
                )
             == kCCMIONoErr
             && intval != 0
            )
            {
                dict.add(opt, intval);
            }
        }

        cellTable_.insert(Id, dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ccm::reader::printInfo() const
{
    Info<< "Mesh Information" << nl
        << "----------------" << nl
        << "boundingBox: " << boundBox(points_) << nl
        << "nPoints: " << nPoints_ << nl
        << "nCells: " << nCells_ << nl
        << "nFaces: " << nFaces_ << nl
        << "nInternalFaces: " << nInternalFaces_ << nl
        << "nBaffles: " << bafInterfaces_.size() << endl;
}


bool Foam::ccm::reader::readGeometry
(
    const scalar scaleFactor
)
{
    detectGeometry();

    if (geometryStatus_ == OKAY)
    {
        // Sanity on the scaling factor

        readMeshTopology(scaleFactor <= VSMALL ? 1 : scaleFactor);
        reorderMesh();

        // assume success if we have points and cells
        if (nCells_ && points_.size())
        {
            geometryStatus_ = READ;
        }
        else
        {
            geometryStatus_ = BAD;
        }
    }

    return (geometryStatus_ == OKAY || geometryStatus_ == READ);
}


bool Foam::ccm::reader::hasGeometry()
{
    detectGeometry();
    return (geometryStatus_ == OKAY || geometryStatus_ == READ);
}


bool Foam::ccm::reader::hasSolution()
{
    detectSolution();
    return (solutionStatus_ == OKAY || solutionStatus_ == READ);
}


void Foam::ccm::reader::writeMesh
(
    const polyMesh& mesh,
    IOstreamOption streamOpt
) const
{
    mesh.removeFiles();

    Info<< "Writing polyMesh" << endl;
    mesh.writeObject(streamOpt, true);
    writeAux(mesh);
}


bool Foam::ccm::reader::remapMeshInfo
(
    const objectRegistry& registry,
    const fileName& remappingDictName
)
{
    dictionary remapDict;

    if (remappingDictName.empty())
    {
        remapDict = IOdictionary
        (
            IOobject
            (
                "remapping",
                "constant",
                registry,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
    }
    else
    {
        // Specified (absolute) name: treat like MUST_READ
        remapDict = dictionary(IFstream(remappingDictName)());
    }

    bool ok = false;

    if (remapDict.empty())
    {
        return false;
    }

    const dictionary* subDictPtr = nullptr;

    // Merge specified cellTable entries together
    if ((subDictPtr = remapDict.findDict("cellTable")) != nullptr)
    {
        cellTable_.combine(*subDictPtr, cellTableId_);
        ok = true;
    }

    // Rename boundaries
    if ((subDictPtr = remapDict.findDict("boundaryRegion")) != nullptr)
    {
        boundaryRegion_.rename(*subDictPtr);
        ok = true;
    }

    return ok;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from reading a file
Foam::ccm::reader::reader(const fileName& file, const reader::options& opts)
:
    base(),
    options_(new options(opts)),
    geometryStatus_(UNKNOWN),
    solutionStatus_(UNKNOWN),
    interfaceDefinitions_(),
    boundaryRegion_(),
    cellTable_(),
    nPoints_(0),
    nInternalFaces_(0),
    nFaces_(0),
    nCells_(0),
    points_(),
    faces_(),
    faceOwner_(),
    faceNeighbour_(),
    bafInterfaces_(),
    domInterfaces_(),
    origFaceId_(),
    origCellId_(),
    cellTableId_(),
    origBndId_(),
    patchSizes_(),
    solutionTable_(),
    fieldTable_(),
    lagrangianTable_()
{
    if (!option().keptSomeRegion())
    {
        FatalErrorInFunction
            << "must retain at least one region type: fluid | porous | solid"
            << exit(FatalError);
    }

    if (!isFile(file, false))
    {
        FatalErrorInFunction
            << "Cannot read file " << file
            << exit(FatalError);
    }

    // Reinitialize error state from the return value
    globalState_->error = CCMIOOpenFile
    (
        nullptr,
        file.c_str(),
        kCCMIORead,
        &(globalState_->root)
    );
    assertNoError("Error opening file for reading");

    detectGeometry();
    detectSolution();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ccm::reader::~reader()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::ccm::reader::options& Foam::ccm::reader::option() const
{
    return *options_;
}


// ************************************************************************* //

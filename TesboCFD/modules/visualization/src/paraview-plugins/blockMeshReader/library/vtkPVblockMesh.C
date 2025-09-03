/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "vtkPVblockMesh.H"
#include "vtkPVblockMeshReader.h"

// OpenFOAM includes
#include "blockMesh.H"
#include "blockMeshTools.H"
#include "foamVersion.H"
#include "Time.H"
#include "OSspecific.H"  // For isDir, cwd
#include "patchZones.H"
#include "StringStream.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkPVblockMesh, 0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    // file-scope
    //- Create a text actor
    static vtkSmartPointer<vtkTextActor> createTextActor
    (
        const std::string& s,
        const Foam::point& pt
    )
    {
        auto txt = vtkSmartPointer<vtkTextActor>::New();

        txt->SetInput(s.c_str());

        // Set text properties
        vtkTextProperty* tprop = txt->GetTextProperty();
        tprop->SetFontFamilyToArial();
        tprop->BoldOn();
        tprop->ShadowOff();
        tprop->SetLineSpacing(1.0);
        tprop->SetFontSize(14);
        tprop->SetColor(1.0, 0.0, 1.0);
        tprop->SetJustificationToCentered();

        txt->GetPositionCoordinate()->SetCoordinateSystemToWorld();
        txt->GetPositionCoordinate()->SetValue(pt.x(), pt.y(), pt.z());

        return txt;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkPVblockMesh::resetCounters()
{
    // Reset mesh part ids and sizes
    rangeBlocks_.reset();
    rangeEdges_.reset();
    rangeCorners_.reset();
}


void Foam::vtkPVblockMesh::updateInfoBlocks
(
    vtkDataArraySelection* select
)
{
    DebugInfo
        << "<beg> updateInfoBlocks"
        << " [meshPtr=" << (meshPtr_ ? "set" : "null") << "]" << endl;

    rangeBlocks_.reset(select->GetNumberOfArrays());

    const blockMesh& blkMesh = *meshPtr_;

    const int nBlocks = blkMesh.size();
    for (int blockI = 0; blockI < nBlocks; ++blockI)
    {
        const blockDescriptor& blockDef = blkMesh[blockI];

        // Display either blockI as a number or with its name
        // (looked up from blockMeshDict)
        OStringStream ostr;
        blockDescriptor::write(ostr, blockI, blkMesh.meshDict());

        // append the (optional) zone name
        if (!blockDef.zoneName().empty())
        {
            ostr << " - " << blockDef.zoneName();
        }

        // Add "blockId" or "blockId - zoneName" to GUI list
        select->AddArray(ostr.str().c_str());
        ++rangeBlocks_;
    }

    DebugInfo << "<end> updateInfoBlocks" << endl;
}


void Foam::vtkPVblockMesh::updateInfoEdges
(
    vtkDataArraySelection* select
)
{
    DebugInfo
        << "<beg> updateInfoEdges"
        << " [meshPtr=" << (meshPtr_ ? "set" : "null") << "]" << endl;

    rangeEdges_.reset(select->GetNumberOfArrays());

    const blockMesh& blkMesh = *meshPtr_;

    for (const blockEdge& blkEdge : blkMesh.edges())
    {
        OStringStream ostr;
        blockVertex::write(ostr, blkEdge.start(), blkMesh.meshDict());
        ostr<< ":";
        blockVertex::write(ostr, blkEdge.end(), blkMesh.meshDict());
        ostr << " - " << blkEdge.type();

        // Add "beg:end - type" to GUI list
        select->AddArray(ostr.str().c_str());
        ++rangeEdges_;
    }

    DebugInfo << "<end> updateInfoEdges" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPVblockMesh::vtkPVblockMesh
(
    const char* const FileName,
    vtkPVblockMeshReader* reader
)
:
    reader_(reader),
    dbPtr_(nullptr),
    meshPtr_(nullptr),
    meshRegion_(polyMesh::defaultRegion),
    meshDir_(polyMesh::meshSubDir),
    rangeBlocks_("block"),
    rangeEdges_("edges"),
    rangeCorners_("corners")
{
    // Reduce some output
    ::Foam::infoDetailLevel = 0;

    DebugInfo << "vtkPVblockMesh - " << FileName << endl;

    // Avoid argList (possible side-effects)
    // - get rootPath/caseName directly from the file
    fileName fullCasePath(fileName(FileName).path());

    if (!isDir(fullCasePath))
    {
        return;
    }
    if (fullCasePath == ".")
    {
        fullCasePath = cwd();
    }

    // OPENFOAM API
    setEnv("FOAM_API", std::to_string(foamVersion::api), true);

    // The name of the executable, unless already present in the environment
    setEnv("FOAM_EXECUTABLE", "paraview", false);

    // Set the case as an environment variable - some BCs might use this
    if (fullCasePath.name().find("processor", 0) == 0)
    {
        const fileName globalCase = fullCasePath.path();

        setEnv("FOAM_CASE", globalCase, true);
        setEnv("FOAM_CASENAME", globalCase.name(), true);
    }
    else
    {
        setEnv("FOAM_CASE", fullCasePath, true);
        setEnv("FOAM_CASENAME", fullCasePath.name(), true);
    }

    // look for 'case{region}.OpenFOAM'
    // could be stringent and insist the prefix match the directory name...
    // Note: cannot use fileName::name() due to the embedded '{}'
    string caseName(fileName(FileName).lessExt());
    string::size_type beg = caseName.find_last_of("/{");
    string::size_type end = caseName.find('}', beg);

    if
    (
        beg != string::npos && caseName[beg] == '{'
     && end != string::npos && end == caseName.size()-1
    )
    {
        meshRegion_ = caseName.substr(beg+1, end-beg-1);

        // some safety
        if (meshRegion_.empty())
        {
            meshRegion_ = polyMesh::defaultRegion;
        }

        if (meshRegion_ != polyMesh::defaultRegion)
        {
            meshDir_ = meshRegion_/polyMesh::meshSubDir;
        }
    }

    DebugInfo
        << "fullCasePath=" << fullCasePath << nl
        << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
        << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << endl;

    // Create time object
    dbPtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            fileName(fullCasePath.path()),
            fileName(fullCasePath.name())
        )
    );

    dbPtr_().functionObjects().off();

    updateInfo();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkPVblockMesh::~vtkPVblockMesh()
{
    DebugInfo << "~vtkPVblockMesh" << endl;

    delete meshPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVblockMesh::updateInfo()
{
    DebugInfo
        << "<beg> updateInfo"
        << " [meshPtr=" << (meshPtr_ ? "set" : "null") << "] " << endl;

    resetCounters();

    vtkDataArraySelection* blockSelection = reader_->GetBlockSelection();
    vtkDataArraySelection* edgeSelection  = reader_->GetCurvedEdgesSelection();

    const bool firstTime = (!blockSelection->GetNumberOfArrays() && !meshPtr_);

    // Preserve the enabled selections if possible
    HashSet<string> enabledParts;
    HashSet<string> enabledEdges;
    if (!firstTime)
    {
        enabledParts = getSelectedArraySet(blockSelection);
        enabledEdges = getSelectedArraySet(edgeSelection);
    }

    // Clear current mesh parts list
    blockSelection->RemoveAllArrays();
    edgeSelection->RemoveAllArrays();

    // need a blockMesh
    updateFoamMesh();

    // Update mesh parts list
    updateInfoBlocks(blockSelection);

    // Update curved edges list
    updateInfoEdges(edgeSelection);

    // Restore the enabled selections
    if (!firstTime)
    {
        setSelectedArrayEntries(blockSelection, enabledParts);
        setSelectedArrayEntries(edgeSelection,  enabledEdges);
    }

    DebugInfo << "<end> updateInfo" << endl;
}


void Foam::vtkPVblockMesh::updateFoamMesh()
{
    DebugInfo << "<beg> updateFoamMesh" << endl;

    // Check to see if the OpenFOAM mesh has been created
    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Creating blockMesh at time=" << dbPtr_().timeName()
                << endl;
        }

        // Set path for the blockMeshDict
        const word dictName("blockMeshDict");
        fileName dictPath(dbPtr_().system()/dictName);

        // Check if dictionary is present in the constant directory
        if
        (
            exists
            (
                dbPtr_().path()/dbPtr_().constant()
               /polyMesh::meshSubDir/dictName
            )
        )
        {
            dictPath = dbPtr_().constant()/polyMesh::meshSubDir/dictName;
        }

        // Store dictionary since is used as database inside blockMesh class
        // for names of vertices and blocks
        IOdictionary* meshDictPtr = new IOdictionary
        (
            IOobject
            (
                dictPath,
                dbPtr_(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );
        meshDictPtr->store();

        meshPtr_ = new blockMesh(*meshDictPtr, meshRegion_);
    }


    DebugInfo << "<end> updateFoamMesh" << endl;
}


void Foam::vtkPVblockMesh::Update
(
    vtkMultiBlockDataSet* output
)
{
    reader_->UpdateProgress(0.1);

    // Update the OpenFOAM mesh
    updateFoamMesh();
    reader_->UpdateProgress(0.5);

    // Convert mesh elemente
    int blockNo = 0;

    convertMeshCorners(output, blockNo);
    convertMeshBlocks(output, blockNo);
    convertMeshEdges(output, blockNo);

    reader_->UpdateProgress(0.8);
}


void Foam::vtkPVblockMesh::UpdateFinalize()
{
    reader_->UpdateProgress(1.0);
}


void Foam::vtkPVblockMesh::renderPatchNames
(
    vtkRenderer* renderer,
    const bool show
)
{
    // always remove old actors first
    forAll(patchTextActors_, actori)
    {
        renderer->RemoveViewProp(patchTextActors_[actori]);
    }
    patchTextActors_.clear();

    // the number of text actors
    label nActors = 0;

    if (show && meshPtr_)
    {
        const blockMesh& blkMesh = *meshPtr_;
        const dictionary& meshDescription = blkMesh.meshDict();

        const pointField cornerPts(blkMesh.vertices(true));

        if (!meshDescription.found("boundary"))
        {
            return;
        }

        // 8 sides per block is plenty
        patchTextActors_.resize(8*blkMesh.size());

        // Collect all variables
        dictionary varDict(meshDescription.subOrEmptyDict("namedVertices"));
        varDict.merge(meshDescription.subOrEmptyDict("namedBlocks"));

        // Read like boundary file
        const PtrList<entry> patchesInfo(meshDescription.lookup("boundary"));

        for (const entry& patchInfo : patchesInfo)
        {
            if (!patchInfo.isDict())
            {
                IOWarningInFunction(meshDescription)
                    << "Entry " << patchInfo << " in boundary section is not a"
                    << " valid dictionary."
                    << endl;
                break;
            }

            const word& patchName = patchInfo.keyword();

            // Read block faces
            faceList patchFaces = blockMeshTools::read<face>
            (
                patchInfo.dict().lookup("faces"),
                varDict
            );

            for (const face& f : patchFaces)
            {
                // Into a list for later removal
                patchTextActors_[nActors++] = createTextActor
                (
                    patchName,
                    f.centre(cornerPts)
                );

                if (nActors == patchTextActors_.size())
                {
                    // hit max allocated space - bail out
                    break;
                }
            }

            if (nActors == patchTextActors_.size())
            {
                // hit max allocated space - bail out
                break;
            }
        }

        patchTextActors_.resize(nActors);
    }

    // Add text to each renderer
    forAll(patchTextActors_, actori)
    {
        renderer->AddViewProp(patchTextActors_[actori]);
    }
}


void Foam::vtkPVblockMesh::renderPointNumbers
(
    vtkRenderer* renderer,
    const bool show
)
{
    // always remove old actors first

    forAll(pointTextActors_, actori)
    {
        renderer->RemoveViewProp(pointTextActors_[actori]);
    }
    pointTextActors_.clear();

    if (show && meshPtr_)
    {
        const blockMesh& blkMesh = *meshPtr_;
        const pointField cornerPts(blkMesh.vertices(true));

        pointTextActors_.resize(cornerPts.size());
        forAll(cornerPts, pointi)
        {
            // Display either pointi as a number or with its name
            // (looked up from blockMeshDict)
            OStringStream os;
            blockVertex::write(os, pointi, blkMesh.meshDict());

            // Into a list for later removal
            pointTextActors_[pointi] = createTextActor
            (
                os.str(),
                cornerPts[pointi]
            );
        }
    }

    // Add text to each renderer
    forAll(pointTextActors_, actori)
    {
        renderer->AddViewProp(pointTextActors_[actori]);
    }
}


void Foam::vtkPVblockMesh::PrintSelf(ostream& os, vtkIndent indent) const
{
#if 0
    os  << indent << "Number of nodes: "
        << (meshPtr_ ? meshPtr_->nPoints() : 0) << "\n";

    os  << indent << "Number of cells: "
        << (meshPtr_ ? meshPtr_->nCells() : 0) << "\n";

    os  << indent << "Number of available time steps: "
        << (dbPtr_ ? dbPtr_->times().size() : 0) << endl;
#endif
}


// ************************************************************************* //

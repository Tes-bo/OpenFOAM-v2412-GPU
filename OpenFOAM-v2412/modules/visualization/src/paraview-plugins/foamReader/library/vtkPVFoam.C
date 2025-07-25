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

#include "vtkPVFoam.H"
#include "vtkPVFoamReader.h"

// OpenFOAM includes
#include "areaFaMesh.H"
#include "faMesh.H"
#include "fvMesh.H"
#include "foamVersion.H"
#include "Time.H"
#include "OSspecific.H"  // For isDir, cwd
#include "patchZones.H"
#include "IOobjectList.H"
#include "collatedFileOperation.H"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"

// Templates (only needed here)
#include "vtkPVFoamUpdateTemplates.C"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkPVFoam, 0);

    // file-scope
    static word updateStateName(polyMesh::readUpdateState state)
    {
        switch (state)
        {
            case polyMesh::UNCHANGED:      return "UNCHANGED";
            case polyMesh::POINTS_MOVED:   return "POINTS_MOVED";
            case polyMesh::TOPO_CHANGE:    return "TOPO_CHANGE";
            case polyMesh::TOPO_PATCH_CHANGE: return "TOPO_PATCH_CHANGE";
        };

        return "UNKNOWN";
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    // file-scope

    //- Create a text actor
    vtkSmartPointer<vtkTextActor> createTextActor
    (
        const std::string& s,
        const Foam::point& pt
    )
    {
        vtkSmartPointer<vtkTextActor> txt =
            vtkSmartPointer<vtkTextActor>::New();

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

void Foam::vtkPVFoam::resetCounters()
{
    // Reset array range information (ids and sizes)
    rangeVolume_.reset();
    rangePatches_.reset();
    rangeClouds_.reset();
    rangeCellZones_.reset();
    rangeFaceZones_.reset();
    rangePointZones_.reset();
    rangeCellSets_.reset();
    rangeFaceSets_.reset();
    rangePointSets_.reset();
}


template<class Container>
bool Foam::vtkPVFoam::addOutputBlock
(
    vtkMultiBlockDataSet* output,
    const HashTable<Container, string>& cache,
    const arrayRange& selector,
    const bool singleDataset
) const
{
    const auto blockNo = output->GetNumberOfBlocks();
    vtkSmartPointer<vtkMultiBlockDataSet> block;
    int datasetNo = 0;

    const List<label> partIds = selector.intersection(selectedPartIds_);

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word shortName = getFoamName(longName);

        auto iter = cache.find(longName);
        if (iter.good() && iter.val().dataset)
        {
            auto dataset = iter.val().dataset;

            if (singleDataset)
            {
                output->SetBlock(blockNo, dataset);
                output->GetMetaData(blockNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    shortName.c_str()
                );

                ++datasetNo;
                break;
            }
            else if (datasetNo == 0)
            {
                block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
                output->SetBlock(blockNo, block);
                output->GetMetaData(blockNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    selector.name()
                );
            }

            block->SetBlock(datasetNo, dataset);
            block->GetMetaData(datasetNo)->Set
            (
                vtkCompositeDataSet::NAME(),
                shortName.c_str()
            );

            ++datasetNo;
        }
    }

    return datasetNo;
}


int Foam::vtkPVFoam::setTime(const std::vector<double>& requestTimes)
{
    if (requestTimes.empty())
    {
        return 0;
    }

    Time& runTime = dbPtr_();

    // Get times list. Flush first to force refresh.
    fileHandler().flush();
    instantList Times = runTime.times();

    int nearestIndex = timeIndex_;
    for (const double& timeValue : requestTimes)
    {
        const int index = Time::findClosestTimeIndex(Times, timeValue);
        if (index >= 0 && index != timeIndex_)
        {
            nearestIndex = index;
            break;
        }
    }

    if (nearestIndex < 0)
    {
        nearestIndex = 0;
    }

    if (debug)
    {
        Info<< "<beg> setTime(";
        unsigned reqi = 0;
        for (const double& timeValue : requestTimes)
        {
            if (reqi) Info<< ", ";
            Info<< timeValue;
            ++reqi;
        }
        Info<< ") - previousIndex = " << timeIndex_
            << ", nearestIndex = " << nearestIndex << nl;
    }

    // See what has changed
    if (timeIndex_ != nearestIndex)
    {
        timeIndex_ = nearestIndex;
        runTime.setTime(Times[nearestIndex], nearestIndex);

        // When mesh changes, so do fields
        meshState_ =
        (
            volMeshPtr_
          ? volMeshPtr_->readUpdate()
          : polyMesh::TOPO_CHANGE
        );

        reader_->UpdateProgress(0.05);

        // this seems to be needed for catching Lagrangian fields
        updateInfo();
    }

    DebugInfo
        << "<end> setTime() - selectedTime="
        << Times[nearestIndex].name() << " index=" << timeIndex_
        << "/" << Times.size()
        << " meshUpdateState=" << updateStateName(meshState_) << nl;

    return nearestIndex;
}


Foam::word Foam::vtkPVFoam::getReaderPartName(const int partId) const
{
    return getFoamName(reader_->GetPartArrayName(partId));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkPVFoam::vtkPVFoam
(
    const char* const vtkFileName,
    vtkPVFoamReader* reader
)
:
    reader_(reader),
    dbPtr_(nullptr),
    volMeshPtr_(nullptr),
    areaMeshPtr_(nullptr),
    meshRegion_(polyMesh::defaultRegion),
    meshDir_(polyMesh::meshSubDir),
    timeIndex_(-1),
    decomposePoly_(false),
    meshState_(polyMesh::TOPO_CHANGE),
    rangeVolume_("volMesh"),
    rangeArea_("areaMesh"),
    rangePatches_("boundary"),
    rangeClouds_("lagrangian"),
    rangeCellZones_("cellZone"),
    rangeFaceZones_("faceZone"),
    rangePointZones_("pointZone"),
    rangeCellSets_("cellSet"),
    rangeFaceSets_("faceSet"),
    rangePointSets_("pointSet")
{
    // Reduce some output
    ::Foam::infoDetailLevel = 0;

    if (debug)
    {
        Info<< "vtkPVFoam - " << vtkFileName << nl;
        printMemory();
    }

    fileName FileName(vtkFileName);

    // Avoid argList (possible side-effects)
    // - get rootPath/caseName directly from the file
    fileName fullCasePath(FileName.path());

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
    if (fullCasePath.name().find("processors", 0) == 0)
    {
        // FileName e.g. "cavity/processors256/processor1.OpenFOAM
        // Remove the processors section so it goes into processorDDD
        // checking below.
        fullCasePath = fullCasePath.path()/fileName(FileName.name()).lessExt();
    }


    if (fullCasePath.name().find("processor", 0) == 0)
    {
        // Give filehandler opportunity to analyse number of processors
        (void)fileHandler().filePath(fullCasePath);

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
    const auto beg = caseName.find_last_of("/{");
    const auto end = caseName.find('}', beg);

    if
    (
        beg != std::string::npos && caseName[beg] == '{'
     && end != std::string::npos && end == caseName.size()-1
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
        << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
        << "region=" << meshRegion_ << nl;

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

Foam::vtkPVFoam::~vtkPVFoam()
{
    DebugInfo << "~vtkPVFoam" << nl;

    delete volMeshPtr_;
    delete areaMeshPtr_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::updateInfo()
{
    DebugInfo
        << "<beg> updateInfo"
        << " [volMeshPtr=" << (volMeshPtr_ ? "set" : "nullptr")
        << "] timeIndex="
        << timeIndex_ << nl;

    resetCounters();

    // Part selection
    {
        vtkDataArraySelection* select = reader_->GetPartSelection();

        // There are two ways to ensure we have the correct list of parts:
        // 1. remove everything and then set particular entries 'on'
        // 2. build a 'char **' list and call SetArraysWithDefault()
        //
        // Nr. 2 has the potential advantage of not touching the modification
        // time of the vtkDataArraySelection, but the qt/paraview proxy
        // layer doesn't care about that anyhow.

        HashSet<string> enabled;
        if (!select->GetNumberOfArrays() && !volMeshPtr_)
        {
            // Fake enable 'internalMesh' on the first call
            enabled = { "internalMesh" };
        }
        else
        {
            // Preserve the enabled selections
            enabled = getSelectedArraySet(select);
        }

        select->RemoveAllArrays();   // Clear existing list

        // Update mesh parts list - add Lagrangian at the bottom
        updateInfoInternalMesh(select);
        updateInfoPatches(select, enabled);
        updateInfoSets(select);
        updateInfoZones(select);
        updateInfoAreaMesh(select);
        updateInfoLagrangian(select);

        setSelectedArrayEntries(select, enabled);  // Adjust/restore selected
    }

    // Volume and area fields - includes save/restore of selected
    updateInfoContinuumFields(reader_->GetVolFieldSelection());

    // Point fields - includes save/restore of selected
    updateInfoPointFields(reader_->GetPointFieldSelection());

    // Lagrangian fields - includes save/restore of selected
    updateInfoLagrangianFields(reader_->GetLagrangianFieldSelection());

    DebugInfo << "<end> updateInfo" << nl;
}


void Foam::vtkPVFoam::Update
(
    vtkMultiBlockDataSet* output,
    vtkMultiBlockDataSet* outputLagrangian
)
{
    if (debug)
    {
        cout<< "<beg> Foam::vtkPVFoam::Update\n";
        output->Print(cout);
        if (outputLagrangian) outputLagrangian->Print(cout);
        printMemory();
    }
    reader_->UpdateProgress(0.1);

    const int caching = reader_->GetMeshCaching();
    const bool oldDecomp = decomposePoly_;
    decomposePoly_ = !reader_->GetUseVTKPolyhedron();

    // Set up mesh parts selection(s)
    // Update cached, saved, unneed values.
    {
        vtkDataArraySelection* selection = reader_->GetPartSelection();
        const int n = selection->GetNumberOfArrays();

        selectedPartIds_.clear();
        HashSet<string> nowActive;

        for (int id=0; id < n; ++id)
        {
            const string str(selection->GetArrayName(id));
            const bool status = selection->GetArraySetting(id);

            if (status)
            {
                selectedPartIds_.set(id, str); // id -> name
                nowActive.set(str);
            }

            if (debug > 1)
            {
                Info<< "  part[" << id << "] = " << status
                    << " : " << str << nl;
            }
        }

        // Dispose of unneeded components
        cachedVtp_.retain(nowActive);
        cachedVtu_.retain(nowActive);

        if
        (
            !caching
         || meshState_ == polyMesh::TOPO_CHANGE
         || meshState_ == polyMesh::TOPO_PATCH_CHANGE
        )
        {
            // Eliminate cached values that would be unreliable
            forAllIters(cachedVtp_, iter)
            {
                iter.val().clearGeom();
                iter.val().clear();
            }
            forAllIters(cachedVtu_, iter)
            {
                iter.val().clearGeom();
                iter.val().clear();
            }
        }
        else if (oldDecomp != decomposePoly_)
        {
            // poly-decompose changed - dispose of cached values
            forAllIters(cachedVtu_, iter)
            {
                iter.val().clearGeom();
                iter.val().clear();
            }
        }
    }

    reader_->UpdateProgress(0.15);

    // Update the OpenFOAM mesh
    {
        if (debug)
        {
            Info<< "<beg> updateFoamMesh" << nl;
            printMemory();
        }

        if (!caching)
        {
            delete volMeshPtr_;
            delete areaMeshPtr_;

            volMeshPtr_ = nullptr;
            areaMeshPtr_ = nullptr;
        }

        // Check to see if the OpenFOAM mesh has been created
        if (!volMeshPtr_)
        {
            DebugInfo
                << "Creating OpenFOAM mesh for region " << meshRegion_
                << " at time=" << dbPtr_().timeName() << nl;

            volMeshPtr_ = new fvMesh
            (
                IOobject
                (
                    meshRegion_,
                    dbPtr_().timeName(),
                    dbPtr_(),
                    IOobject::MUST_READ
                )
            );

            meshState_ = polyMesh::TOPO_CHANGE; // New mesh
        }
        else
        {
            DebugInfo << "Using existing OpenFOAM mesh" << nl;
        }

        if (rangeArea_.intersects(selectedPartIds_))
        {
            if (!areaMeshPtr_)
            {
                areaMeshPtr_ = new faMesh(*volMeshPtr_);
            }
        }
        else
        {
            delete areaMeshPtr_;

            areaMeshPtr_ = nullptr;
        }

        if (debug)
        {
            Info<< "<end> updateFoamMesh" << nl;
            printMemory();
        }
    }

    reader_->UpdateProgress(0.4);

    convertMeshVolume();
    convertMeshPatches();
    reader_->UpdateProgress(0.6);

    if (reader_->GetIncludeZones())
    {
        convertMeshCellZones();
        convertMeshFaceZones();
        convertMeshPointZones();
        reader_->UpdateProgress(0.65);
    }

    if (reader_->GetIncludeSets())
    {
        convertMeshCellSets();
        convertMeshFaceSets();
        convertMeshPointSets();
        reader_->UpdateProgress(0.7);
    }

    convertMeshArea();

    convertMeshLagrangian();

    reader_->UpdateProgress(0.8);

    // Update fields
    convertVolFields();
    convertPointFields();
    convertAreaFields();

    convertLagrangianFields();

    // Assemble multiblock output
    addOutputBlock(output, cachedVtu_, rangeVolume_, true); // One dataset
    addOutputBlock(output, cachedVtp_, rangePatches_);
    addOutputBlock(output, cachedVtu_, rangeCellZones_);
    addOutputBlock(output, cachedVtp_, rangeFaceZones_);
    addOutputBlock(output, cachedVtp_, rangePointZones_);
    addOutputBlock(output, cachedVtu_, rangeCellSets_);
    addOutputBlock(output, cachedVtp_, rangeFaceSets_);
    addOutputBlock(output, cachedVtp_, rangePointSets_);
    addOutputBlock(output, cachedVtp_, rangeArea_);
    addOutputBlock
    (
        (outputLagrangian ? outputLagrangian : output),
        cachedVtp_,
        rangeClouds_
    );

    DebugInfo << "done reader part" << nl << nl;
    reader_->UpdateProgress(0.95);

    meshState_ = polyMesh::UNCHANGED;

    if (caching & 2)
    {
        // Suppress caching of Lagrangian since it normally always changes.
        cachedVtp_.filterKeys
        (
            [](const word& k){ return k.starts_with("lagrangian/"); },
            true // prune
        );
    }
    else
    {
        cachedVtp_.clear();
        cachedVtu_.clear();
    }
}


void Foam::vtkPVFoam::UpdateFinalize()
{
    if (!reader_->GetMeshCaching())
    {
        delete volMeshPtr_;
        delete areaMeshPtr_;

        volMeshPtr_ = nullptr;
        areaMeshPtr_ = nullptr;
    }

    reader_->UpdateProgress(1.0);
}


std::vector<double> Foam::vtkPVFoam::findTimes(const bool skipZero) const
{
    std::vector<double> times;

    if (dbPtr_)
    {
        const Time& runTime = dbPtr_();
        // Get times list. Flush first to force refresh.
        fileHandler().flush();
        instantList timeLst = runTime.times();

        // find the first time for which this mesh appears to exist
        label begIndex = timeLst.size();
        forAll(timeLst, timei)
        {
            if
            (
                IOobject
                (
                    "points",
                    timeLst[timei].name(),
                    meshDir_,
                    runTime
                ).typeHeaderOk<pointIOField>(false, false)
            )
            {
                begIndex = timei;
                break;
            }
        }

        label nTimes = timeLst.size() - begIndex;

        // skip "constant" time whenever possible
        if (begIndex == 0 && nTimes > 1)
        {
            if (timeLst[begIndex].name() == runTime.constant())
            {
                ++begIndex;
                --nTimes;
            }
        }

        // skip "0/" time if requested and possible
        if (skipZero && nTimes > 1 && timeLst[begIndex].name() == "0")
        {
            ++begIndex;
            --nTimes;
        }

        times.reserve(nTimes);
        while (nTimes-- > 0)
        {
            times.push_back(timeLst[begIndex++].value());
        }
    }
    else
    {
        if (debug)
        {
            cout<< "no valid dbPtr:\n";
        }
    }

    return times;
}


void Foam::vtkPVFoam::renderPatchNames
(
    vtkRenderer* renderer,
    const bool show
)
{
    // Always remove old actors first

    for (auto& actor : patchTextActors_)
    {
        renderer->RemoveViewProp(actor);
    }
    patchTextActors_.clear();

    if (show && volMeshPtr_)
    {
        // Get the display patches, strip off any prefix/suffix
        wordHashSet selectedPatches
        (
            getSelected(reader_->GetPartSelection(), rangePatches_)
        );

        if (selectedPatches.empty())
        {
            return;
        }

        const polyBoundaryMesh& pbm = volMeshPtr_->boundaryMesh();

        // Find the total number of zones
        // Each zone will take the patch name
        // Number of zones per patch ... zero zones should be skipped
        labelList nZones(pbm.size(), Zero);

        // Per global zone number the average face centre position
        List<DynamicList<point>> zoneCentre(pbm.size());


        // Loop through all patches to determine zones, and centre of each zone
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            // Only include the patch if it is selected
            if (!selectedPatches.contains(pp.name()))
            {
                continue;
            }

            const labelListList& edgeFaces = pp.edgeFaces();
            const vectorField& n = pp.faceNormals();

            boolList featEdge(pp.nEdges(), false);

            forAll(edgeFaces, edgei)
            {
                const labelList& eFaces = edgeFaces[edgei];

                if (eFaces.size() == 1)
                {
                    // Note: could also do ones with > 2 faces but this gives
                    // too many zones for baffles
                    featEdge[edgei] = true;
                }
                else if (mag(n[eFaces[0]] & n[eFaces[1]]) < 0.5)
                {
                    featEdge[edgei] = true;
                }
            }

            // Do topological analysis of patch, find disconnected regions
            patchZones pZones(pp, featEdge);

            nZones[patchi] = pZones.nZones();

            labelList zoneNFaces(pZones.nZones(), Zero);

            // Create storage for additional zone centres
            forAll(zoneNFaces, zonei)
            {
                zoneCentre[patchi].append(Zero);
            }

            // Do averaging per individual zone
            forAll(pp, facei)
            {
                const label zonei = pZones[facei];
                zoneCentre[patchi][zonei] += pp[facei].centre(pp.points());
                zoneNFaces[zonei]++;
            }

            forAll(zoneCentre[patchi], zonei)
            {
                zoneCentre[patchi][zonei] /= zoneNFaces[zonei];
            }
        }

        // Count number of zones we're actually going to display.
        // This is truncated to a max per patch

        const label MAXPATCHZONES = 20;

        label displayZoneI = 0;

        forAll(pbm, patchi)
        {
            displayZoneI += min(MAXPATCHZONES, nZones[patchi]);
        }

        DebugInfo
            << "displayed zone centres = " << displayZoneI << nl
            << "zones per patch = " << nZones << nl;

        // Set the size of the patch labels to max number of zones
        patchTextActors_.resize(displayZoneI);

        DebugInfo << "constructing patch labels" << nl;

        // Actor index
        displayZoneI = 0;

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            // Only selected patches will have a non-zero number of zones
            const label nDisplayZones = min(MAXPATCHZONES, nZones[patchi]);
            label increment = 1;
            if (nZones[patchi] >= MAXPATCHZONES)
            {
                increment = nZones[patchi]/MAXPATCHZONES;
            }

            label globalZoneI = 0;
            for (label i = 0; i < nDisplayZones; ++i, globalZoneI += increment)
            {
                DebugInfo
                    << "patch name = " << pp.name() << nl
                    << "anchor = " << zoneCentre[patchi][globalZoneI] << nl
                    << "globalZoneI = " << globalZoneI << nl;

                // Into a list for later removal
                patchTextActors_[displayZoneI++] = createTextActor
                (
                    pp.name(),
                    zoneCentre[patchi][globalZoneI]
                );
            }
        }

        // Resize the patch names list to the actual number of patch names added
        patchTextActors_.resize(displayZoneI);
    }

    // Add text to each renderer
    for (auto& actor : patchTextActors_)
    {
        renderer->AddViewProp(actor);
    }
}


void Foam::vtkPVFoam::PrintSelf(ostream& os, vtkIndent indent) const
{
    os  << indent << "Number of nodes: "
        << (volMeshPtr_ ? volMeshPtr_->nPoints() : 0) << "\n";

    os  << indent << "Number of cells: "
        << (volMeshPtr_ ? volMeshPtr_->nCells() : 0) << "\n";

    os  << indent << "Number of available time steps: "
        << (dbPtr_ ? dbPtr_->times().size() : 0) << "\n";

    os  << indent << "mesh region: " << meshRegion_ << "\n";
}


void Foam::vtkPVFoam::printInfo() const
{
    std::cout
        << "Region: " << meshRegion_ << "\n"
        << "nPoints: " << (volMeshPtr_ ? volMeshPtr_->nPoints() : 0) << "\n"
        << "nCells: "  << (volMeshPtr_ ? volMeshPtr_->nCells() : 0) << "\n"
        << "nTimes: "
        << (dbPtr_ ? dbPtr_->times().size() : 0) << "\n";

    std::vector<double> times = this->findTimes(reader_->GetSkipZeroTime());

    std::cout<<"  " << times.size() << "(";
    for (const double& val : times)
    {
        std::cout<< ' ' << val;
    }
    std::cout << " )" << std::endl;
}


// ************************************************************************* //

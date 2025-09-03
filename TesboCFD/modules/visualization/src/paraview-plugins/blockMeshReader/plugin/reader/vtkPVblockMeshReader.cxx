/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "vtkPVblockMeshReader.h"

#include "pqApplicationCore.h"
#include "pqRenderView.h"
#include "pqServerManagerModel.h"

// VTK includes
#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkSMRenderViewProxy.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"

// OpenFOAM includes
#include "vtkPVblockMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

vtkStandardNewMacro(vtkPVblockMeshReader);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vtkPVblockMeshReader::vtkPVblockMeshReader()
{
    Debug = 0;
    vtkDebugMacro(<<"Constructor");

    SetNumberOfInputPorts(0);

    FileName = nullptr;
    backend_ = nullptr;

    ShowPatchNames   = false;
    ShowPointNumbers = true;

    BlockSelection = vtkDataArraySelection::New();
    CurvedEdgesSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback
    (
        &vtkPVblockMeshReader::SelectionModifiedCallback
    );
    SelectionObserver->SetClientData(this);

    BlockSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );

    CurvedEdgesSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vtkPVblockMeshReader::~vtkPVblockMeshReader()
{
    vtkDebugMacro(<<"Destructor");

    if (backend_)
    {
        // Remove text actors
        updatePatchNamesView(false);
        updatePointNumbersView(false);

        delete backend_;
        backend_ = nullptr;
    }

    if (FileName)
    {
        delete[] FileName;
    }

    BlockSelection->RemoveAllObservers();
    CurvedEdgesSelection->RemoveAllObservers();

    SelectionObserver->Delete();
    BlockSelection->Delete();
    CurvedEdgesSelection->Delete();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

int vtkPVblockMeshReader::RequestInformation
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestInformation");

    if (!FileName)
    {
        vtkErrorMacro("FileName has to be specified!");
        return 0;
    }

    if (Foam::vtkPVblockMesh::debug)
    {
        cout<<"REQUEST_INFORMATION\n";
        outputVector->GetInformationObject(0)->Print(cout);
    }

    if (backend_)
    {
        backend_->updateInfo();
    }
    else
    {
        backend_ = new Foam::vtkPVblockMesh(FileName, this);
    }

    return 1;
}


int vtkPVblockMeshReader::RequestData
(
    vtkInformation* vtkNotUsed(request),
    vtkInformationVector** vtkNotUsed(inputVector),
    vtkInformationVector* outputVector
)
{
    vtkDebugMacro(<<"RequestData");

    if (!FileName)
    {
        vtkErrorMacro("FileName must be specified!");
        return 0;
    }
    if (!backend_)
    {
        // Catch some previous error
        vtkErrorMacro("Reader failed - perhaps no blockMesh?");
        return 0;
    }

    if (Foam::vtkPVblockMesh::debug)
    {
        cout<<"RequestData:\n";
        outputVector->GetInformationObject(0)->Print(cout);
    }

    auto* output = vtkMultiBlockDataSet::SafeDownCast
    (
        outputVector->GetInformationObject(0)->Get
        (
            vtkMultiBlockDataSet::DATA_OBJECT()
        )
    );

    backend_->Update(output);

    updatePatchNamesView(ShowPatchNames);
    updatePointNumbersView(ShowPointNumbers);

    backend_->UpdateFinalize();

    return 1;
}


void vtkPVblockMeshReader::Refresh()
{
    // Delete the current blockMesh to force re-read and update
    if (backend_)
    {
        // Remove text actors
        updatePatchNamesView(false);
        updatePointNumbersView(false);

        delete backend_;
        backend_ = nullptr;
    }

    this->Modified();
}


void vtkPVblockMeshReader::SetShowPatchNames(bool val)
{
    if (ShowPatchNames != val)
    {
        ShowPatchNames = val;
        updatePatchNamesView(ShowPatchNames);
    }
}


void vtkPVblockMeshReader::SetShowPointNumbers(bool val)
{
    if (ShowPointNumbers != val)
    {
        ShowPointNumbers = val;
        updatePointNumbersView(ShowPointNumbers);
    }
}


void vtkPVblockMeshReader::updatePatchNamesView(const bool show)
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    // Need to check this, since our destructor calls this
    if (!appCore)
    {
        return;
    }

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();
    if (!smModel || !backend_)
    {
        return;
    }

    // Get all the pqRenderView instances
    for (auto view : smModel->findItems<pqRenderView*>())
    {
        backend_->renderPatchNames
        (
            view->getRenderViewProxy()->GetRenderer(),
            show
        );
    }

    // Use refresh here?
}


void vtkPVblockMeshReader::updatePointNumbersView(const bool show)
{
    pqApplicationCore* appCore = pqApplicationCore::instance();

    // Need to check this, since our destructor calls this
    if (!appCore)
    {
        return;
    }

    // Server manager model for querying items in the server manager
    pqServerManagerModel* smModel = appCore->getServerManagerModel();
    if (!smModel || !backend_)
    {
        return;
    }

    // Get all the pqRenderView instances
    for (auto view : smModel->findItems<pqRenderView*>())
    {
        backend_->renderPointNumbers
        (
            view->getRenderViewProxy()->GetRenderer(),
            show
        );
    }

    // Use refresh here?
}


void vtkPVblockMeshReader::PrintSelf(ostream& os, vtkIndent indent)
{
    vtkDebugMacro(<<"PrintSelf");

    this->Superclass::PrintSelf(os,indent);
    os  << indent << "File name: "
        << (this->FileName ? this->FileName : "(none)") << "\n";

    backend_->PrintSelf(os, indent);
}


// ----------------------------------------------------------------------
// Block selection list control

vtkDataArraySelection* vtkPVblockMeshReader::GetBlockSelection()
{
    return BlockSelection;
}

int vtkPVblockMeshReader::GetNumberOfBlockArrays()
{
    return BlockSelection->GetNumberOfArrays();
}

const char* vtkPVblockMeshReader::GetBlockArrayName(int index)
{
    return BlockSelection->GetArrayName(index);
}

int vtkPVblockMeshReader::GetBlockArrayStatus(const char* name)
{
    return BlockSelection->ArrayIsEnabled(name);
}

void vtkPVblockMeshReader::SetBlockArrayStatus
(
    const char* name,
    int status
)
{
    if (status)
    {
        BlockSelection->EnableArray(name);
    }
    else
    {
        BlockSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------
// CurvedEdges selection list control

vtkDataArraySelection* vtkPVblockMeshReader::GetCurvedEdgesSelection()
{
    return CurvedEdgesSelection;
}

int vtkPVblockMeshReader::GetNumberOfCurvedEdgesArrays()
{
    return CurvedEdgesSelection->GetNumberOfArrays();
}

const char* vtkPVblockMeshReader::GetCurvedEdgesArrayName(int index)
{
    return CurvedEdgesSelection->GetArrayName(index);
}

int vtkPVblockMeshReader::GetCurvedEdgesArrayStatus(const char* name)
{
    return CurvedEdgesSelection->ArrayIsEnabled(name);
}

void vtkPVblockMeshReader::SetCurvedEdgesArrayStatus
(
    const char* name,
    int status
)
{
    if (status)
    {
        CurvedEdgesSelection->EnableArray(name);
    }
    else
    {
        CurvedEdgesSelection->DisableArray(name);
    }
}


// ----------------------------------------------------------------------

void vtkPVblockMeshReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkPVblockMeshReader*>(clientdata)->Modified();
}


int vtkPVblockMeshReader::FillOutputPortInformation
(
    int port,
    vtkInformation* info
)
{
    if (port == 0)
    {
        return this->Superclass::FillOutputPortInformation(port, info);
    }
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
}


// ************************************************************************* //

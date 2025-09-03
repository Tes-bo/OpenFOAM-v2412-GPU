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

Class
    vtkPVblockMeshReader

Description
    reads a dataset in OpenFOAM bockMesh format

    vtkPVblockMeshReader creates an multiblock dataset.
    It uses the OpenFOAM infrastructure (blockMesh).

SourceFiles
    vtkPVblockMeshReader.cxx

\*---------------------------------------------------------------------------*/

#ifndef vtkPVblockMeshReader_h
#define vtkPVblockMeshReader_h

// VTK includes
#include "vtkMultiBlockDataSetAlgorithm.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

// VTK forward declarations
class vtkDataArraySelection;
class vtkCallbackCommand;

namespace Foam
{
    class vtkPVblockMesh;
}

/*---------------------------------------------------------------------------*\
                   Class vtkPVblockMeshReader Declaration
\*---------------------------------------------------------------------------*/

class vtkPVblockMeshReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeMacro(vtkPVblockMeshReader, vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream&, vtkIndent) override;

    static vtkPVblockMeshReader* New();

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // Display patch names
    virtual void SetShowPatchNames(bool);
    vtkGetMacro(ShowPatchNames, bool);

    // Description:
    // Display corner point labels
    virtual void SetShowPointNumbers(bool);
    vtkGetMacro(ShowPointNumbers, bool);

    // Description:
    // Refresh blockMesh from changes to blockMeshDict
    virtual void Refresh();

    // Description:
    // Blocks selection list control
    vtkDataArraySelection* GetBlockSelection();
    int  GetNumberOfBlockArrays();
    int  GetBlockArrayStatus(const char* name);
    void SetBlockArrayStatus(const char* name, int status);
    const char* GetBlockArrayName(int index);

    // Description:
    // CurvedEdges selection list control
    vtkDataArraySelection* GetCurvedEdgesSelection();
    int  GetNumberOfCurvedEdgesArrays();
    int  GetCurvedEdgesArrayStatus(const char* name);
    void SetCurvedEdgesArrayStatus(const char* name, int status);
    const char* GetCurvedEdgesArrayName(int index);

    // Description:
    // Callback registered with the SelectionObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );


protected:

    //- Default construct
    vtkPVblockMeshReader();

    //- Destructor
    virtual ~vtkPVblockMeshReader();

    //- Return information about mesh, times, etc without loading anything
    virtual int RequestInformation
    (
        vtkInformation* unusedRequest,
        vtkInformationVector** unusedInputVector,
        vtkInformationVector* outputVector
    ) override;

    //- Get the mesh for a particular time
    virtual int RequestData
    (
        vtkInformation* unusedRequest,
        vtkInformationVector** unusedInputVector,
        vtkInformationVector* outputVector
    ) override;

    //- Fill in additional port information
    virtual int FillOutputPortInformation(int, vtkInformation*) override;

    // The observer to modify this object when array selections are modified
    vtkCallbackCommand* SelectionObserver;

    char* FileName;


private:

    //- No copy construct
    vtkPVblockMeshReader(const vtkPVblockMeshReader&) = delete;

    //- No copy assignment
    void operator=(const vtkPVblockMeshReader&) = delete;

    //- Add/remove patch names to/from the view
    void updatePatchNamesView(const bool show);

    //- Add/remove point numbers to/from the view
    void updatePointNumbersView(const bool show);


    //- Show Patch Names
    bool ShowPatchNames;

    //- Show Point Numbers
    bool ShowPointNumbers;

    vtkDataArraySelection* BlockSelection;
    vtkDataArraySelection* CurvedEdgesSelection;

    //- Backend portion of the reader
    Foam::vtkPVblockMesh* backend_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

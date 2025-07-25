/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::vtkPVblockMesh

Description
    The backend for the vtkPVblockMeshReader reader module -
    providing a paraview reader interface for OpenFOAM blockMesh.

    The block reader module can assist when creating a blockMeshDict
    for use with the blockMesh utility. As well as blocks, it can be
    used to visualize edges,corners and patch names.

    There is no native VTK equivalent for this functionality.

SourceFiles
    vtkPVblockMesh.C
    vtkPVblockMeshConvert.C

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPVblockMesh_H
#define vtkPVblockMesh_H

#include "foamPvCore.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkPoints;
class vtkPVblockMeshReader;
class vtkRenderer;
class vtkTextActor;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkUnstructuredGrid;
class vtkIndent;

template<class T> class vtkSmartPointer;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations (OpenFOAM)
class Time;
class blockMesh;
template<class Type> class List;

/*---------------------------------------------------------------------------*\
                     Class vtkPVblockMesh Declaration
\*---------------------------------------------------------------------------*/

class vtkPVblockMesh
:
    private foamPvCore
{
    // Private Data

        //- Access to the controlling vtkPVblockMeshReader
        vtkPVblockMeshReader* reader_;

        //- OpenFOAM time control
        autoPtr<Time> dbPtr_;

        //- OpenFOAM mesh
        blockMesh* meshPtr_;

        //- The mesh region
        word meshRegion_;

        //- The mesh directory for the region
        fileName meshDir_;

        //- First instance and size of bleckMesh blocks
        //  used to index into blockStatus_
        arrayRange rangeBlocks_;

        //- First instance and size of CurvedEdges (only partially used)
        arrayRange rangeEdges_;

        //- First instance and size of block corners (only partially used)
        arrayRange rangeCorners_;

        //- List of patch names for rendering to window
        List<vtkSmartPointer<vtkTextActor>> patchTextActors_;

        //- List of point numbers for rendering to window
        List<vtkSmartPointer<vtkTextActor>> pointTextActors_;


    // Private Member Functions

        //- Reset data counters
        void resetCounters();

        //- OpenFOAM mesh
        void updateFoamMesh();

        //- Internal block info
        void updateInfoBlocks(vtkDataArraySelection* select);

        //- Block curved edges info
        void updateInfoEdges(vtkDataArraySelection* select);

        //- Mesh blocks
        void convertMeshBlocks(vtkMultiBlockDataSet*, int& blockNo);

        //- Mesh curved edges
        void convertMeshEdges(vtkMultiBlockDataSet*, int& blockNo);

        //- Mesh corners
        void convertMeshCorners(vtkMultiBlockDataSet*, int& blockNo);


        //- No copy construct
        vtkPVblockMesh(const vtkPVblockMesh&) = delete;

        //- No copy assignment
        void operator=(const vtkPVblockMesh&) = delete;


public:

    //- Declare type-name (with debug switch)
    ClassName("vtkPVblockMesh");


    // Constructors

        //- Construct from components
        vtkPVblockMesh
        (
            const char* const FileName,
            vtkPVblockMeshReader* reader
        );


    //- Destructor
    ~vtkPVblockMesh();


    // Member Functions

        //- Update
        void updateInfo();

        void Update(vtkMultiBlockDataSet* output);

        //- Final part of Update(), after any last minute rendering.
        void UpdateFinalize();

        //- Add/remove patch names to/from the view
        void renderPatchNames(vtkRenderer* renderer, const bool show);

        //- Add/remove point numbers to/from the view
        void renderPointNumbers(vtkRenderer* renderer, const bool show);


    // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

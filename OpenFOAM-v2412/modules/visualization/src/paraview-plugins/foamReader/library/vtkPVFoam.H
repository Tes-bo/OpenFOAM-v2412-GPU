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
    Foam::vtkPVFoam

Description
    The backend for the vtkPVFoamReader reader module -
    providing a paraview reader interface for OpenFOAM meshes and fields.

    Similar, and sometimes better, functionality may be provided by the
    native VTK OpenFOAM reader. OpenCFD has recently (2017) been working
    on improving the native VTK OpenFOAM reader for the benefit of everyone.

    In some areas the reader module lacks compared to the native reader
    (notably the ability to work on decomosed datasets), but provides
    additional handling of sets,zones,groups. Some features have also since
    been adapted to the native reader. Additionally, the reader module
    provides a useful platform for testing new ideas.

Note
    The reader module allows two levels of caching. The OpenFOAM fvMesh
    can be cached in memory, for faster loading of fields. Additionally,
    the translated VTK geometries are held in a local cache. The cached
    VTK geometries should incur no additional overhead since they use
    the VTK reference counting for their storage management.

SourceFiles
    vtkPVFoam.C
    vtkPVFoamFields.C
    vtkPVFoamMesh.C
    vtkPVFoamMeshLagrangian.C
    vtkPVFoamMeshVolume.C
    vtkPVFoamTemplates.C
    vtkPVFoamUpdateInfo.C
    vtkPVFoamFieldTemplates.C
    vtkPVFoamUpdateTemplates.C

    // Needed by VTK:
    vtkDataArrayTemplateImplicit.txx

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoam_H
#define vtkPVFoam_H

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "wordList.H"
#include "primitivePatch.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "foamPvCore.H"
#include "foamVtkVtuAdaptor.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkCellArray;
class vtkDataArraySelection;
class vtkDataSet;
class vtkFloatArray;
class vtkIndent;
class vtkPVFoamReader;
class vtkRenderer;
class vtkTextActor;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations (OpenFOAM)
class Time;
class faMesh;
class fvMesh;
class IOobjectList;
class polyPatch;
class fvMeshSubset;

template<class Type> class IOField;
template<class Type> class Field;
template<class Type> class List;

/*---------------------------------------------------------------------------*\
                        Class vtkPVFoam Declaration
\*---------------------------------------------------------------------------*/

class vtkPVFoam
:
    private foamPvCore
{
    // Convenience typedefs
    typedef PrimitivePatchInterpolation<primitivePatch> patchInterpolator;

        //- Bookkeeping for vtkPolyData
        struct foamVtpData
        :
            public vtk::Caching<vtkPolyData>,
            public foamVtkMeshMaps
        {};

        //- Bookkeeping for vtkUnstructuredGrid
        struct foamVtuData
        :
            public vtk::vtuAdaptor
        {
            //- Subsetted mesh as vtkUnstructuredGrid
            vtkSmartPointer<vtkUnstructuredGrid> subset
            (
                const fvMeshSubset& subsetter,
                bool decompPoly = false
            );
        };


    // Private Data

        //- Access to the controlling vtkPVFoamReader
        vtkPVFoamReader* reader_;

        //- OpenFOAM time control
        autoPtr<Time> dbPtr_;

        //- OpenFOAM finite volume mesh
        fvMesh* volMeshPtr_;

        //- OpenFOAM finite area mesh
        faMesh* areaMeshPtr_;

        //- The mesh region
        word meshRegion_;

        //- The mesh directory for the region
        fileName meshDir_;

        //- The time index
        int timeIndex_;

        //- Previous/current decomposition request
        bool decomposePoly_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- The index of selected parts mapped to their names
        Map<string> selectedPartIds_;

        //- Any information for 2D (VTP) geometries
        HashTable<foamVtpData, string> cachedVtp_;

        //- Cell maps and other information for 3D (VTU) geometries
        HashTable<foamVtuData, string> cachedVtu_;

        //- First instance and size of various mesh parts
        //  used to index into selectedPartIds and thus indirectly into
        //  cachedVtp, cachedVtu
        arrayRange rangeVolume_;
        arrayRange rangeArea_;
        arrayRange rangePatches_;
        arrayRange rangeClouds_;
        arrayRange rangeCellZones_;
        arrayRange rangeFaceZones_;
        arrayRange rangePointZones_;
        arrayRange rangeCellSets_;
        arrayRange rangeFaceSets_;
        arrayRange rangePointSets_;

        //- List of patch names for rendering to window
        List<vtkSmartPointer<vtkTextActor>> patchTextActors_;


    // Private Member Functions

        template<class Container>
        bool addOutputBlock
        (
            vtkMultiBlockDataSet* output,
            const HashTable<Container, string>& cache,
            const arrayRange& selector,
            const bool singleDataset = false
        ) const;


        //- Reset data counters
        void resetCounters();


    // Update information helper functions

        //- Internal mesh info
        void updateInfoInternalMesh(vtkDataArraySelection* select);

        //- Finite area mesh info
        void updateInfoAreaMesh(vtkDataArraySelection* select);

        //- Lagrangian info
        void updateInfoLagrangian(vtkDataArraySelection* select);

        //- Patch info, modifies enabledEntries
        void updateInfoPatches
        (
            vtkDataArraySelection* select,
            HashSet<string>& enabledEntries
        );

        //- Set info
        void updateInfoSets(vtkDataArraySelection* select);

        //- Zone info
        void updateInfoZones(vtkDataArraySelection* select);


        //- Get non-empty zone names for zoneType from file
        wordList getZoneNames(const word& zoneType) const;

        //- Get names of on non-empty zones from the mesh info
        template<class ZoneType>
        static wordList getZoneNames
        (
            const ZoneMesh<ZoneType, polyMesh>& zmesh
        );

        //- Field info
        template<template<class> class patchType, class meshType>
        void updateInfoFields
        (
            vtkDataArraySelection* select,
            const IOobjectList& objects
        );

        //- Volume/Area field info
        void updateInfoContinuumFields(vtkDataArraySelection* select);

        //- Point field info
        void updateInfoPointFields(vtkDataArraySelection* select);

        //- Lagrangian field info
        void updateInfoLagrangianFields(vtkDataArraySelection* select);


    // Mesh conversion functions

        //- Convert internalMesh
        void convertMeshVolume();

        //- Convert areaMesh
        void convertMeshArea();

        //- Convert Lagrangian points
        void convertMeshLagrangian();

        //- Convert mesh patches.
        //  The additionalIds (cached data) contain the patch Ids.
        //  There will be several for groups, but only one for regular patches.
        void convertMeshPatches();

        //- Convert cell zones
        void convertMeshCellZones();

        //- Convert cell sets
        void convertMeshCellSets();

        //- Convert face zones
        void convertMeshFaceZones();

        //- Convert face sets
        //  The cellMap (cached data) contains the face-labels.
        void convertMeshFaceSets();

        //- Convert point zones
        //  The pointMap (cached data) contains the point-labels.
        void convertMeshPointZones();

        //- Convert point sets
        //  The pointMap (cached data) contains the point-labels.
        void convertMeshPointSets();


    // Add mesh functions

        //- Lagrangian positions as vtkPolyData
        vtkSmartPointer<vtkPolyData> lagrangianVTKMesh
        (
            const polyMesh& mesh,
            const word& cloudName
        ) const;


    // Field conversion functions

        //- Face set/zone field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertFaceFieldToVTK
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const labelUList& faceLabels
        ) const;


        //- Convert finite volume fields
        void convertVolFields();

        //- Convert finite area fields
        void convertAreaFields();

        //- Convert point fields
        void convertPointFields();

        //- Convert Lagrangian fields
        void convertLagrangianFields();


    // Convert OpenFOAM fields

        //- Volume field - all types
        template<class Type>
        void convertVolField
        (
            const PtrList<patchInterpolator>& patchInterpList,
            const GeometricField<Type, fvPatchField, volMesh>& fld
        );

        //- Volume fields - all types
        template<class Type>
        void convertVolFields
        (
            const fvMesh& mesh,
            const PtrList<patchInterpolator>& patchInterpList,
            const IOobjectList& objects
        );

        //- Volume internal fields (DimensionedField)- all types
        template<class Type>
        void convertDimFields
        (
            const fvMesh& mesh,
            const PtrList<patchInterpolator>& patchInterpList,
            const IOobjectList& objects
        );

        //- Area fields - all types
        template<class Type>
        void convertAreaFields
        (
            const faMesh& mesh,
            const IOobjectList& objects
        );

        //- Volume field - all selected parts
        template<class Type>
        void convertVolFieldBlock
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr,
            const arrayRange& range
        );

        //- Lagrangian fields - all types
        template<class Type>
        void convertLagrangianFields
        (
            const IOobjectList& objects,
            vtkPolyData* vtkmesh
        );

        //- Point fields - all types
        template<class Type>
        void convertPointFields
        (
            const pointMesh& pMesh,
            const IOobjectList& objectst
        );

        //- Point field - all selected parts
        template<class Type>
        void convertPointFieldBlock
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            const arrayRange& range
        );

        //- Point field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertPointField
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            const GeometricField<Type, fvPatchField, volMesh>& vfld,
            const foamVtuData& vtuData
        );


    // GUI selection helper functions

        //- Get the first word from the reader 'parts' selection
        word getReaderPartName(const int partId) const;


    // Constructors

        //- No copy construct
        vtkPVFoam(const vtkPVFoam&) = delete;

        //- No copy assignment
        void operator=(const vtkPVFoam&) = delete;


public:

    //- Declare type-name (with debug switch)
    ClassName("vtkPVFoam");


    // Constructors

        //- Construct from components
        vtkPVFoam
        (
            const char* const FileName,
            vtkPVFoamReader* reader
        );


    //- Destructor
    ~vtkPVFoam();


    // Member Functions

        //- Update
        void updateInfo();

        void Update
        (
            vtkMultiBlockDataSet* output,
            vtkMultiBlockDataSet* outputLagrangian
        );

        //- Final part of Update(), after any last minute rendering.
        void UpdateFinalize();

        //- Add/remove patch names to/from the view
        void renderPatchNames(vtkRenderer* renderer, const bool show);

        //- Return a list of selected times.
        //  Use STL container since these values are used by the plugin
        std::vector<double> findTimes(const bool skipZero = false) const;

        //- Set the runTime to the first plausible request time,
        //  returns the timeIndex
        //  sets to "constant" on error
        int setTime(const std::vector<double>& requestTimes);


        //- The current time index
        int timeIndex() const
        {
           return timeIndex_;
        }


    // Access

        //- Debug information
        void PrintSelf(ostream&, vtkIndent) const;

        void printInfo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

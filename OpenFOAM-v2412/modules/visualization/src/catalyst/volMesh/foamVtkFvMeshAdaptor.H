/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::vtk::fvMeshAdaptor

Description
    The backend for the catalystFvMesh function object.

    The output is a multi-block dataset with two blocks corresponding
    to the internal (volume) mesh (block 0) and the boundary (block 1),
    which are further divided into sub-blocks for each patch.

    The lowest data blocks are multi-piece datasets (UnstructuredGrid
    and PolyData for internal mesh and boundary, respectively) with
    each piece corresponding to its MPI rank.

Note
    It is possible to restrict the output of the blocks (eg, patches only)
    by specifying the appropriate "channels" prior to conversion.

SourceFiles
    foamVtkFvMeshAdaptor.C
    foamVtkFvMeshAdaptorGeom.C
    foamVtkFvMeshAdaptorGeomVtu.C
    foamVtkFvMeshAdaptorFields.C
    foamVtkFvMeshAdaptorFieldTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamvtk_fvMeshAdaptor_H
#define foamvtk_fvMeshAdaptor_H

#include "className.H"
#include "wordList.H"
#include "Enum.H"
#include "primitivePatch.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"

#include "foamVtkVtuAdaptor.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkCellArray;
class vtkDataSet;
class vtkFloatArray;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
template<class Type> class Field;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::fvMeshAdaptor Declaration
\*---------------------------------------------------------------------------*/

class fvMeshAdaptor
{
public:

    // Public Data Members

        //- The types of output blocks. Defined as a bitmask.
        enum channelType : unsigned
        {
            NONE = 0,           //!< No output (placeholder).
            INTERNAL = 0x1,     //!< Internal mesh.
            BOUNDARY = 0x2,     //!< Boundary patches.
            ALL = 0x3,          //!< Both internal mesh and boundary patches.
            DEFAULT = ALL       //!< Default is ALL
        };

        static const Enum<channelType> channelNames;

        //- Name for internal mesh ("internal")
        static const word& internalName()
        {
            return channelNames[channelType::INTERNAL];
        }

        //- Name for the boundary ("boundary")
        static const word& boundaryName()
        {
            return channelNames[channelType::BOUNDARY];
        }


private:

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
        {};


    // Private Data

        //- OpenFOAM mesh
        const fvMesh& mesh_;

        //- In non-empty, restrict to these selected patches only
        wordRes patchPatterns_;

        //- Cached values for selected patches. *Always* in sorted order.
        labelList patchIds_;

        //- Selected output channels
        unsigned channels_;

        //- Interpolate volume to point fields (default: true)
        //  Generates PointData in VTK
        bool interpFields_;

        //- Extrapolate the internal field to the boundaries (default: false)
        bool extrapPatches_;

        //- Previous/current decomposition request (default: false)
        bool decomposePoly_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- Any information for 2D (VTP) geometries
        HashTable<foamVtpData, string> cachedVtp_;

        //- Cell maps and other information for 3D (VTU) geometries
        HashTable<foamVtuData, string> cachedVtu_;


    // Mesh Conversion

        //- Define patch ids
        void definePatchIds();

        //- Convert internal
        void convertGeometryInternal();

        //- Convert patches
        //  Patches may use additionalIds (cached data) with the patch Ids.
        //  There will be several for groups, but only one for regular patches.
        void convertGeometryBoundary();

        //- Ghost cells - internal mesh
        void applyGhostingInternal(const labelUList& types);

        //- Ghost cells - patches
        void applyGhostingBoundary(const labelUList& types);

        //- Ghost cells to affect the visibility of the geometry.
        //  Currently only used for overset.
        void applyGhosting();


    // Field Conversion

        //- Convert specified volume fields
        void convertVolFields(const wordRes& selectFields);

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
            const PtrList<patchInterpolator>& patchInterpList,
            const wordRes& selectFields
        );

        //- Volume internal fields (DimensionedField)- all types
        template<class Type>
        void convertDimFields
        (
            const PtrList<patchInterpolator>& patchInterpList,
            const wordRes& selectFields
        );

        //- Volume field - all selected parts
        template<class Type>
        void convertVolFieldInternal
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr
        );

        //- Point field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertPointField
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            const GeometricField<Type, fvPatchField, volMesh>& vfld,
            const foamVtuData& vtuData
        );


        //- Update geometry and fields
        void updateContent(const wordRes& selectFields);


    // Constructors

        //- No copy construct
        fvMeshAdaptor(const fvMeshAdaptor&) = delete;

        //- No copy assignment
        void operator=(const fvMeshAdaptor&) = delete;


public:

    //- Declare type-name (with debug switch)
    ClassName("vtk::fvMeshAdaptor");


    // Constructors

        //- Construct from components
        explicit fvMeshAdaptor
        (
            const fvMesh& mesh,
            const channelType channelsOpt = channelType::DEFAULT,
            const wordRes& patchSelection = wordRes()
        );


    //- Destructor
    ~fvMeshAdaptor() = default;


    // Member Functions

        //- Define the output channels by name
        void setChannels(const wordList& chanNames);

        //- Define the output channels by enum
        void setChannels(enum channelType chanIds);

        //- Define the output channels by value
        void setChannels(unsigned chanIds);

        //- Define polyhedral decomposition treatment
        void setDecompose(const bool on);


        //- Return the selected output channel ids
        label channels() const;

        //- True if INTERNAL channel is being used
        bool usingInternal() const;

        //- True if BOUNDARY channel is being used
        bool usingBoundary() const;

        //- Selected (non-processor) patch ids, when the BOUNDARY channel
        //- is being used. Empty otherwise.
        const labelList& patchIds() const;

        //- Return the names of known (supported) fields
        wordHashSet knownFields(const wordRes& selectFields) const;

        void updateState(polyMesh::readUpdateState state);

        //- The output is a pair (INTERNAL/BOUNDARY channels) of vtk meshes
        //- with point/cell data.
        vtkSmartPointer<vtkMultiBlockDataSet> output
        (
            const wordRes& selectFields
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

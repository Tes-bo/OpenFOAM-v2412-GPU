/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Class
    Foam::functionObjects::vtkWrite

Group
    grpUtilitiesFunctionObjects

Description
    Writes fields in VTK (xml or legacy) format.
    Writes cell-values or point-interpolated values for volFields.

    Example of function object specification:
    \verbatim
    vtkWrite1
    {
        type            vtkWrite;
        libs            (utilityFunctionObjects);
        writeControl    writeTime;
        writeInterval   1;
        format          binary;
        legacy          false;
        decompose       false;
        ...
        fields          (U p);
        // excludeFields  ("force.*");

        selection
        {
            box
            {
                action  use;
                source  box;
                box     (-0.1 -0.01 -0.1) (0.1 0.30 0.1);
            }
            dome
            {
                action  add;
                source  sphere;
                origin  (-0.1 -0.01 -0.1);
                radius  0.25;
            }
            centre
            {
                action  subtract;
                source  sphere;
                origin  (-0.1 -0.01 -0.1);
                radius  0.1;
            }
            blob
            {
                action  add;
                source  surface;
                surface triSurfaceMesh;
                name    blob.stl;
            }
        }
    }
    \endverbatim

    \heading Basic Usage
    \table
        Property    | Description                           | Required | Default
        type        | Type name: vtkWrite                   | yes |
        fields      | Select fields to output (wordRe list) | yes |
        excludeFields | Exclude fields from output (wordRe list) | no  |
        boundary    | Convert boundary fields               | no  | true
        internal    | Convert internal fields               | no  | true
        single      | Combine patches into a single boundary | no | false
        interpolate | Interpolate for point values          | no  | false
    \endtable

    \heading Output Options
    \table
        Property    | Description                           | Required | Default
        format      | ascii or binary format                | no  | binary
        legacy      | Legacy VTK output                     | no  | false
        precision   | Write precision in ascii         | no | same as IOstream
        directory   | The output directory name     | no | postProcessing/NAME
        width       | Padding width for file name           | no  | 8
        decompose   | Decompose polyhedral cells            | no  | false
        writeIds    | Write cell,patch,proc id fields       | no  | false
    \endtable

    \heading Output Selection
    \table
        Property    | Description                           | Required | Default
        region      | Name for a single region              | no  | region0
        regions     | List of regions (wordRe list)         | no  |
        patches     | Limit to listed patches (wordRe list) | no  |
        excludePatches | Exclude specified patches          | no  |
        selection   | Cell selection (topoSet actions)      | no  | empty dict
    \endtable

Note
    The region of interest is defined by the selection dictionary
    as a set of actions (use,add,subtract,subset,invert).
    Omitting the selection dictionary is the same as specifying the
    conversion of all cells (in the selected regions).
    Omitting the patches entry is the same as specifying the conversion of all
    patches.

See also
    Foam::functionObjects::ensightWrite
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl
    Foam::cellBitSet::select

SourceFiles
    vtkWrite.C
    vtkWriteImpl.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_vtkWrite_H
#define Foam_functionObjects_vtkWrite_H

#include "timeFunctionObject.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"
#include "foamVtkSeriesWriter.H"
#include "fvMeshSubsetProxy.H"
#include "searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class vtkWrite Declaration
\*---------------------------------------------------------------------------*/

class vtkWrite
:
    public functionObjects::timeFunctionObject
{
    // Private Data

        //- The output directory
        fileName outputDir_;

        //- The printf format for zero-padding names
        string printf_;

        //- VTK output options
        vtk::outputOptions writeOpts_;

        //- Verbose output
        bool verbose_;

        //- Convert internal mesh
        bool doInternal_;

        //- Convert boundary mesh
        bool doBoundary_;

        //- Combine patches into a single file
        bool oneBoundary_;

        //- Interpolate cell to point values
        bool interpolate_;

        //- Decompose polyhedra
        bool decompose_;

        //- Write cell ids field
        bool writeIds_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- Requested names of regions to process
        wordRes selectRegions_;

        //- Requested names of patches to process
        wordRes selectPatches_;

        //- Selection of patches to block
        wordRes blockPatches_;

        //- Requested selection of fields to process
        wordRes selectFields_;

        //- Selection of fields to block
        wordRes blockFields_;

        //- Dictionary of volume selections
        dictionary selection_;

        //- Pointers to the requested mesh regions (sorted)
        UPtrList<const fvMesh> meshes_;

        //- Subsetting for meshes.
        //  Access index according to sorted mesh names.
        PtrList<fvMeshSubset> meshSubsets_;

        //- Storage for VTU cells, sizing.
        //  Access index according to sorted mesh names.
        PtrList<vtk::vtuCells> vtuMappings_;

        //- VTK file series
        HashTable<vtk::seriesWriter, fileName> series_;


    // Private Member Functions

        //- Update mesh subset according to zones, geometry, bounds
        bool updateSubset(fvMeshSubset& subsetter) const;

        //- Get patchIds selected in list
        labelList getSelectedPatches(const polyBoundaryMesh& patches) const;

        //- Read information for selections
        bool readSelection(const dictionary& dict);

        //- Update meshes, subMeshes etc.
        bool update();


    // Write

        //- Write all volume fields
        label writeAllVolFields
        (
            autoPtr<vtk::internalWriter>& internalWriter,
            UPtrList<vtk::patchWriter>& patchWriters,
            const fvMeshSubset& proxy,
            const wordHashSet& candidateNames
        ) const;

        //- Write all volume fields with point interpolation
        label writeAllVolFields
        (
            autoPtr<vtk::internalWriter>& internalWriter,
            const autoPtr<volPointInterpolation>& pInterp,
            UPtrList<vtk::patchWriter>& patchWriters,
            const UPtrList
            <
                PrimitivePatchInterpolation<primitivePatch>
            >& patchInterps,
            const fvMeshSubset& proxy,
            const wordHashSet& candidateNames
        ) const;

        //- Write selected GeoField fields.
        template<class GeoField>
        label writeVolFieldsImpl
        (
            autoPtr<vtk::internalWriter>& internalWriter,
            UPtrList<vtk::patchWriter>& patchWriters,
            const fvMeshSubset& proxy,
            const wordHashSet& candidateNames
        ) const;

        //- Write selected GeoField fields with point interpolation
        template<class GeoField>
        label writeVolFieldsImpl
        (
            autoPtr<vtk::internalWriter>& internalWriter,
            const autoPtr<volPointInterpolation>& pInterp,
            UPtrList<vtk::patchWriter>& patchWriters,
            const UPtrList
            <
                PrimitivePatchInterpolation<primitivePatch>
            >& patchInterps,
            const fvMeshSubset& proxy,
            const wordHashSet& candidateNames
        ) const;


        //- No copy construct
        vtkWrite(const vtkWrite&) = delete;

        //- No copy assignment
        void operator=(const vtkWrite&) = delete;


public:

    //- Runtime type information
    TypeName("vtkWrite");


    // Constructors

        //- Construct from Time and dictionary
        vtkWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~vtkWrite() = default;


    // Member Functions

        //- Read the vtkWrite specification
        virtual bool read(const dictionary& dict);

        //- Execute - does nothing
        virtual bool execute();

        //- Write fields
        virtual bool write();

        //- On end - cleanup internal allocations
        virtual bool end();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

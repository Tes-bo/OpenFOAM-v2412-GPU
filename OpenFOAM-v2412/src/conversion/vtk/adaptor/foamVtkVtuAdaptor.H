/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::vtk::vtuAdaptor

Description
    A low-level backend adaptor for converting OpenFOAM volume meshes/fields
    to VTK internal representation.
    The output is a single piece vtkUnstructuredGrid dataset.
    Multiblock composition is done by the caller.

Note
    This was originally part of the catalystFvMesh function object backend,
    which was prototyped in the ParaView reader module.
    This adaptor unifies many common elements of both.

SourceFiles
    foamVtkVtuAdaptorI.H
    foamVtkVtuAdaptorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtuAdaptor_H
#define Foam_vtk_vtuAdaptor_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "foamVtkTools.H"
#include "foamVtkMeshMaps.H"
#include "foamVtuSizing.H"

#include "vtkUnstructuredGrid.h"
#include "vtkMultiBlockDataSet.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkDataSet;
class vtkFloatArray;

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                       Class vtk::vtuAdaptor Declaration
\*---------------------------------------------------------------------------*/

//- Bookkeeping for vtkUnstructuredGrid
class vtuAdaptor
:
    public vtk::Caching<vtkUnstructuredGrid>,
    public foamVtkMeshMaps
{
public:

    // Member Functions

    // Mesh Conversion

        //- The vtk points for the mesh (and decomposition)
        inline vtkSmartPointer<vtkPoints> points
        (
            const fvMesh& mesh
        ) const;

        //- The vtk points for the mesh (and decomposition)
        //- using the provided pointMap
        inline vtkSmartPointer<vtkPoints> points
        (
            const fvMesh& mesh,
            const labelUList& pointMap
        ) const;

        //- Internal mesh as vtkUnstructuredGrid.
        //  Calling this automatically populates the sizing and vtk::Caching
        inline vtkSmartPointer<vtkUnstructuredGrid> internal
        (
            const fvMesh& mesh,
            const bool decompPoly = false
        );


    // Field Conversion

        //- Convert internal volume field (CellData)
        template<class Type>
        static vtkSmartPointer<vtkFloatArray>
        convertField
        (
            const DimensionedField<Type, volMesh>& fld,
            const vtuAdaptor& vtuData
        );

        //- Convert internal volume field (CellData)
        template<class Type>
        static vtkSmartPointer<vtkFloatArray>
        convertField
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const vtuAdaptor& vtuData
        );


        //- Convert internal volume field (CellData)
        template<class Type>
        vtkSmartPointer<vtkFloatArray>
        convertField
        (
            const DimensionedField<Type, volMesh>& fld
        ) const;

        //- Convert internal volume field (CellData)
        template<class Type>
        vtkSmartPointer<vtkFloatArray>
        convertField
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkVtuAdaptorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkVtuAdaptorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

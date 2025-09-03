/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::vtk::faMeshAdaptor

Description
    The backend for the catalystFaMesh function object.

    The output is a multi-piece PolyData dataset corresponding to the
    finite area mesh. Each piece corresponds to its MPI rank.

SourceFiles
    foamVtkFaMeshAdaptor.C
    foamVtkFaMeshAdaptorGeom.C
    foamVtkFaMeshAdaptorFields.C
    foamVtkFaMeshAdaptorFieldTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkFaMeshAdaptor_H
#define foamVtkFaMeshAdaptor_H

#include "className.H"
#include "polyMesh.H"
#include "areaFieldsFwd.H"
#include "foamVtkTools.H"
#include "foamVtkMeshMaps.H"

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiPieceDataSet.h>

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkCellArray;
class vtkDataSet;
class vtkFloatArray;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class faMesh;
template<class Type> class Field;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::faMeshAdaptor Declaration
\*---------------------------------------------------------------------------*/

class faMeshAdaptor
{
public:

    // Public Data Members

        //- Name for internal mesh ("internal")
        static const word internalName;


private:

    // Private Classes

        //- Bookkeeping for vtkPolyData
        struct foamVtpData
        :
            public vtk::Caching<vtkPolyData>,
            public foamVtkMeshMaps
        {};


    // Private Data

        //- OpenFOAM finite-area mesh
        const faMesh& mesh_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- Any information for 2D (VTP) geometries
        HashTable<foamVtpData, string> cachedVtp_;


    // Mesh Conversion

        //- Convert internal
        void convertGeometryInternal();


    // Field Conversion

        //- Convert specified area fields
        void convertAreaFields(const wordRes& selectFields);

        //- Area field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertAreaFieldToVTK
        (
            const GeometricField<Type, faPatchField, areaMesh>& fld,
            const foamVtpData& vtpData
        ) const;


        //- Update geometry and fields
        void updateContent(const wordRes& selectFields);


    // Convert OpenFOAM fields

        //- Area field - all types
        template<class Type>
        void convertAreaField
        (
            const GeometricField<Type, faPatchField, areaMesh>& fld
        );

        //- Area fields - all types
        template<class Type>
        void convertAreaFields(const wordRes& selectFields);


    // Constructors

        //- No copy construct
        faMeshAdaptor(const faMeshAdaptor&) = delete;

        //- No copy assignment
        void operator=(const faMeshAdaptor&) = delete;


public:

    //- Declare type-name (with debug switch)
    ClassName("vtk::faMeshAdaptor");


    // Constructors

        //- Construct for mesh region
        explicit faMeshAdaptor(const faMesh& mesh);


    //- Destructor
    ~faMeshAdaptor() = default;


    // Member Functions

        //- Return the names of known (supported) fields
        wordHashSet knownFields(const wordRes& selectFields) const;

        void updateState(polyMesh::readUpdateState state);

        //- The output with cell data.
        vtkSmartPointer<vtkMultiPieceDataSet> output
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

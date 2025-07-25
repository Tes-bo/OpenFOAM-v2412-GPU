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

Namespace
    Foam::vtk::Tools

Description
    A collection of static methods to assist converting OpenFOAM data
    structures into VTK internal data structures.

    Remapping of the symmTensor order is required in input or output
    directions. OpenFOAM uses (XX, XY, XZ, YY, YZ, ZZ) order,
    VTK uses (XX, YY, ZZ, XY, YZ, XZ) order.

Note
    The class is implemented as headers-only.

SourceFiles
    foamVtkToolsI.H
    foamVtkToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtkTools_H
#define Foam_vtk_vtkTools_H

#include "faceList.H"
#include "pointField.H"
#include "symmTensor.H"
#include "MinMax.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkIdTypeArray.h"
#include "vtkSmartPointer.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkDataSet;
class vtkCellData;
class vtkPointData;

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                        Class vtk::Caching Declaration
\*---------------------------------------------------------------------------*/

//- Bookkeeping for internal caching.
//  Retain an original copy of the geometry as well as a shallow copy
//  with the output fields.
//  The original copy is reused for different timesteps etc.
template<class DataType>
struct Caching
{
    typedef DataType dataType;

    //- The geometry, without any cell/point data
    vtkSmartPointer<dataType> vtkgeom;

    //- The shallow-copy of geometry, plus additional data
    vtkSmartPointer<dataType> dataset;

    //- Number of points associated with the geometry
    inline uint64_t nPoints() const
    {
        return vtkgeom ? vtkgeom->GetNumberOfPoints() : 0;
    }

    //- Clear geometry and dataset
    void clearGeom()
    {
        vtkgeom = nullptr;
        dataset = nullptr;
    }

    //- Return a shallow copy of vtkgeom for manipulation
    vtkSmartPointer<dataType> getCopy() const
    {
        auto copy = vtkSmartPointer<dataType>::New();
        if (vtkgeom)
        {
            copy->ShallowCopy(vtkgeom);
        }
        return copy;
    }

    //- Make a shallow copy of vtkgeom into dataset
    void reuse()
    {
        dataset = vtkSmartPointer<dataType>::New();
        if (vtkgeom)
        {
            dataset->ShallowCopy(vtkgeom);
        }
    }

    //- Set the geometry and make a shallow copy to dataset
    void set(vtkSmartPointer<dataType> geom)
    {
        vtkgeom = geom;
        reuse();
    }

    //- Report basic information to output
    void PrintSelf(std::ostream& os) const
    {
        os << "geom" << nl;
        if (vtkgeom)
        {
            vtkgeom->PrintSelf(std::cout, vtkIndent(2));
        }
        else
        {
            os << "nullptr";
        }
        os << nl;

        os << "copy" << nl;
        if (dataset)
        {
            dataset->PrintSelf(std::cout, vtkIndent(2));
        }
        else
        {
            os << "nullptr";
        }
        os << nl;
    }
};


/*---------------------------------------------------------------------------*\
                             Namespace vtk::Tools
\*---------------------------------------------------------------------------*/

namespace Tools
{
    //- Wrap vtkUnsignedCharArray as a UList
    inline UList<uint8_t> asUList
    (
        vtkUnsignedCharArray* array,
        const label size
    );

    //- Wrap vtkIdTypeArray as a UList
    inline UList<vtkIdType> asUList
    (
        vtkIdTypeArray* array,
        const label size
    );

    //- Return a list of points as vtkPoints
    inline vtkSmartPointer<vtkPoints> Points
    (
        const UList<point>& pts
    );

    //- Return an indirect list of points as vtkPoints
    inline vtkSmartPointer<vtkPoints> Points
    (
        const UList<point>& pts,
        const labelUList& addr
    );

    //- Convert a list of faces (or triFaces) to vtk polygon cells
    template<class Face>
    vtkSmartPointer<vtkCellArray> Faces(const UList<Face>& faces);

    //- Return vtkPolyData of vertices for each point
    inline vtkSmartPointer<vtkPolyData> Vertices
    (
        const UList<point>& pts
    );

    //- Return vtkPolyData of vertices for each point
    inline vtkSmartPointer<vtkPolyData> Vertices
    (
        const UList<point>& pts,
        const labelUList& addr
    );

    //- Min/Max of scalar, or mag() of non-scalars. Includes nullptr check.
    inline scalarMinMax rangeOf(vtkDataArray* data);


    //- Convert OpenFOAM patch to vtkPolyData
    struct Patch
    {
        //- Return local patch points as vtkPoints
        template<class PatchType>
        static vtkSmartPointer<vtkPoints> points(const PatchType& p);

        //- Convert patch faces to vtk polygon cells
        template<class PatchType>
        static vtkSmartPointer<vtkCellArray> faces(const PatchType& p);

        //- Convert patch points/faces to vtkPolyData
        template<class PatchType>
        static vtkSmartPointer<vtkPolyData> mesh(const PatchType& p);

        //- Convert patch face normals to vtkFloatArray
        template<class PatchType>
        static vtkSmartPointer<vtkFloatArray> faceNormals(const PatchType& p);

        //- Return patch face centres as vtkPoints
        template<class PatchType>
        static vtkSmartPointer<vtkPoints> faceCentres(const PatchType& p);

        //- Convert points/faces component to vtkPolyData
        template<class Face>
        static vtkSmartPointer<vtkPolyData> mesh
        (
            const UList<point>& pts,
            const UList<Face>& fcs
        );
    };


    //- Remapping for some OpenFOAM data types (eg, symmTensor)
    //  \param data[in,out] The data to be remapped in-place
    template<class Type>
    inline void remapTuple(float data[]) {}

    //- Template specialization for symmTensor ordering
    template<>
    inline void remapTuple<symmTensor>(float data[]);

    //- Remapping for some OpenFOAM data types (eg, symmTensor)
    //  \param data[in,out] The data to be remapped in-place
    template<class Type>
    inline void remapTuple(double data[]) {}

    //- Template specialization for symmTensor ordering
    template<>
    inline void remapTuple<symmTensor>(double data[]);

    //- Copy/transcribe OpenFOAM data types to VTK format
    //  This allows a change of data type (float vs double) as well as
    //  addressing any swapping issues (eg, symmTensor)
    //
    //  \param output[out] The output scratch space. Must be long enough
    //     to hold the result.
    //  \param val[in] The input data to copy/transcribe
    template<class Type>
    inline void foamToVtkTuple(float output[], const Type& val);

    //- Copy/transcribe OpenFOAM data types to VTK format
    //  This allows a change of data type (float vs double) as well as
    //  addressing any swapping issues (eg, symmTensor)
    //
    //  \param output[out] The output scratch space. Must be long enough
    //     to hold the result.
    //  \param val[in] The input data to copy/transcribe
    template<class Type>
    inline void foamToVtkTuple(double output[], const Type& val);


    // Field Conversion Functions

    //- Copy list to pre-allocated vtk array.
    //  \return number of input items copied
    template<class Type>
    label transcribeFloatData
    (
        vtkFloatArray* array,
        const UList<Type>& input,
        vtkIdType start = 0         //!< The write offset into output array
    );

    //- Create named field initialized to zero
    template<class Type>
    vtkSmartPointer<vtkFloatArray> zeroField
    (
        const word& name,
        const label size
    );

    //- Convert field data to a vtkFloatArray
    template<class Type>
    vtkSmartPointer<vtkFloatArray> convertFieldToVTK
    (
        const word& name,
        const UList<Type>& fld
    );

    //- An identity list of VTK_VERTEX
    inline vtkSmartPointer<vtkCellArray> identityVertices
    (
        const label size
    );

} // End namespace Tools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Specializations

//- Template specialization for symmTensor ordering
template<>
void Foam::vtk::Tools::remapTuple<Foam::symmTensor>(float data[])
{
    std::swap(data[1], data[3]);    // swap XY <-> YY
    std::swap(data[2], data[5]);    // swap XZ <-> ZZ
}


//- Template specialization for symmTensor ordering
template<>
void Foam::vtk::Tools::remapTuple<Foam::symmTensor>(double data[])
{
    std::swap(data[1], data[3]);    // swap XY <-> YY
    std::swap(data[2], data[5]);    // swap XZ <-> ZZ
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkToolsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamVtkToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

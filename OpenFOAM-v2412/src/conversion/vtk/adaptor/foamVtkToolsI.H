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

\*---------------------------------------------------------------------------*/

#include <numeric>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Foam::UList<uint8_t> Foam::vtk::Tools::asUList
(
    vtkUnsignedCharArray* array,
    const label size
)
{
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(size);

    return UList<uint8_t>(array->WritePointer(0, size), size);
}


inline Foam::UList<vtkIdType> Foam::vtk::Tools::asUList
(
    vtkIdTypeArray* array,
    const label size
)
{
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(size);

    return UList<vtkIdType>(array->WritePointer(0, size), size);
}


inline vtkSmartPointer<vtkPoints>
Foam::vtk::Tools::Points(const UList<point>& pts)
{
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    vtkpoints->SetNumberOfPoints(pts.size());

    vtkIdType pointId = 0;
    for (const point& pt : pts)
    {
        vtkpoints->SetPoint(pointId++, pt.v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkPoints>
Foam::vtk::Tools::Points(const UList<point>& pts, const labelUList& addr)
{
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    vtkpoints->SetNumberOfPoints(addr.size());

    vtkIdType pointId = 0;
    for (const label pointi : addr)
    {
        vtkpoints->SetPoint(pointId++, pts[pointi].v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkPolyData>
Foam::vtk::Tools::Vertices(const UList<point>& pts)
{
    auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

    vtkmesh->SetPoints(Tools::Points(pts));
    vtkmesh->SetVerts(Tools::identityVertices(pts.size()));

    return vtkmesh;
}


inline vtkSmartPointer<vtkPolyData>
Foam::vtk::Tools::Vertices(const UList<point>& pts, const labelUList& addr)
{
    auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

    vtkmesh->SetPoints(Tools::Points(pts, addr));
    vtkmesh->SetVerts(Tools::identityVertices(addr.size()));

    return vtkmesh;
}


inline Foam::scalarMinMax Foam::vtk::Tools::rangeOf(vtkDataArray* data)
{
    double range[2]{GREAT, -GREAT};

    if (data)
    {
        if (data->GetNumberOfComponents() == 1)
        {
            data->GetRange(range);
        }
        else
        {
            // Mag
            data->GetRange(range, -1);
        }
    }

    return scalarMinMax(range[0], range[1]);
}


inline vtkSmartPointer<vtkCellArray> Foam::vtk::Tools::identityVertices
(
    const label size
)
{
    auto cells = vtkSmartPointer<vtkCellArray>::New();

    #ifdef VTK_CELL_ARRAY_V2

    // Offsets
    // [0, n1, n1+n2, n1+n2+n3... ]

    auto offsets = vtkSmartPointer<vtkIdTypeArray>::New();

    {
        const vtkIdType nOffsets(size+1);

        offsets->SetNumberOfTuples(nOffsets);

        vtkIdType* iter = offsets->WritePointer(0, nOffsets);

        std::iota(iter, (iter + nOffsets), vtkIdType(0));
    }


    auto connect = vtkSmartPointer<vtkIdTypeArray>::New();

    // Connectivity
    {
        const vtkIdType nConnect(size);

        connect->SetNumberOfTuples(nConnect);

        vtkIdType* iter = connect->WritePointer(0, nConnect);

        std::iota(iter, (iter + nConnect), vtkIdType(0));
    }

    // Move into a vtkCellArray

    cells->SetData(offsets, connect);

    #else

    // In VTK-8.2.0 and older,
    // sizes are interwoven (prefixed) in the connectivity

    // Connectivity size, with prefixed size information

    // Cell connectivity for vertex
    // [size, ids.., size, ids...]  -> therefore  [1, id, 1, id, ...]

    const vtkIdType nElem(size);
    const vtkIdType nConnect(2*size);

    {
        cells->GetData()->SetNumberOfTuples(nConnect);

        vtkIdType* iter = cells->WritePointer(nElem, nConnect);

        // Fill in the connectivity array, with prefixed size information

        for (vtkIdType id = 0; id < nElem; ++id)
        {
            *(iter++) = 1;
            *(iter++) = id;
        }
    }

    #endif

    return cells;
};


template<class Type>
inline void Foam::vtk::Tools::foamToVtkTuple
(
    float output[],
    const Type& val
)
{
    for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        output[cmpt] = component(val, cmpt);
    }
    remapTuple<Type>(output);
}


template<class Type>
inline void Foam::vtk::Tools::foamToVtkTuple
(
    double output[],
    const Type& val
)
{
    for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        output[cmpt] = component(val, cmpt);
    }
    remapTuple<Type>(output);
}


// ************************************************************************* //

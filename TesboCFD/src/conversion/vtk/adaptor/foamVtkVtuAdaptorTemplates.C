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

\*---------------------------------------------------------------------------*/

// VTK includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const DimensionedField<Type, volMesh>& fld,
    const vtuAdaptor& vtuData
)
{
    const int nComp(pTraits<Type>::nComponents);
    const labelUList& cellMap = vtuData.cellMap();

    auto data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetName(fld.name().c_str());
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(cellMap.size());

    // DebugInfo
    //     << "Convert field: " << fld.name()
    //     << " size=" << cellMap.size()
    //     << " (" << fld.size() << " + "
    //     << (cellMap.size() - fld.size())
    //     << ") nComp=" << nComp << endl;


    float scratch[pTraits<Type>::nComponents];

    vtkIdType celli = 0;
    for (const label meshCelli : cellMap)
    {
        vtk::Tools::foamToVtkTuple(scratch, fld[meshCelli]);
        data->SetTuple(celli++, scratch);
    }

    return data;
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const vtuAdaptor& vtuData
)
{
    return convertField<Type>(fld.internalField(), vtuData);
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const DimensionedField<Type, volMesh>& fld
) const
{
    return convertField<Type>(fld, *this);
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    return convertField<Type>(fld, *this);
}


// ************************************************************************* //

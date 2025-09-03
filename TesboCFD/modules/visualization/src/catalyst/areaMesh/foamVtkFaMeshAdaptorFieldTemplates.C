/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#ifndef foamVtkFaMeshAdaptorFieldTemplates_C
#define foamVtkFaMeshAdaptorFieldTemplates_C

// OpenFOAM includes
#include "error.H"

// vtk includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::faMeshAdaptor::convertAreaFields
(
    const wordRes& selectFields
)
{
    const auto& obr = mesh_.thisDb();

    typedef GeometricField<Type, faPatchField, areaMesh> fieldType;

    for (const auto& fld : obr.csorted<fieldType>(selectFields))
    {
        convertAreaField(fld);
    }
}


template<class Type>
void Foam::vtk::faMeshAdaptor::convertAreaField
(
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    const auto& longName = internalName;

    auto iter = cachedVtp_.find(longName);
    if (!iter.good() || !iter.val().dataset)
    {
        // Should not happen, but for safety require a vtk geometry
        return;
    }
    foamVtpData& vtpData = iter.val();
    auto dataset = vtpData.dataset;

    vtkSmartPointer<vtkFloatArray> cdata = convertAreaFieldToVTK
    (
        fld,
        vtpData
    );
    dataset->GetCellData()->AddArray(cdata);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// low-level conversions
//

template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::faMeshAdaptor::convertAreaFieldToVTK
(
    const GeometricField<Type, faPatchField, areaMesh>& fld,
    const foamVtpData& vtpData
) const
{
    // The vtpData is not used for anything
    return vtk::Tools::convertFieldToVTK(fld.name(), fld);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

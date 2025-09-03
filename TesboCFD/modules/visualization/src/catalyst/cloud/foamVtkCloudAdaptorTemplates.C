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

#ifndef foamVtkCloudAdaptorTemplates_C
#define foamVtkCloudAdaptorTemplates_C

#include "foamVtkTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::label Foam::vtk::cloudAdaptor::convertLagrangianFields
(
    vtkPolyData* vtkmesh,
    const objectRegistry& obr
)
{
    typedef IOField<Type> fieldType;

    label nFields = 0;

    for (const auto& fld : obr.csorted<fieldType>())
    {
        vtkSmartPointer<vtkFloatArray> data = vtk::Tools::convertFieldToVTK
        (
            fld.name(),
            fld
        );

        // Provide identical data as cell and as point data
        vtkmesh->GetCellData()->AddArray(data);
        vtkmesh->GetPointData()->AddArray(data);

        ++nFields;
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "vtkPVFoam.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkSmartPointer<vtkUnstructuredGrid> Foam::vtkPVFoam::foamVtuData::subset
(
    const fvMeshSubset& subsetter,
    bool decompPoly
)
{
    auto vtkmesh = internal(subsetter.subMesh(), decompPoly);

    // Convert cellMap, addPointCellLabels to global cell ids
    renumberCells(subsetter.cellMap());

    // Copy pointMap as well, otherwise pointFields fail
    pointMap() = subsetter.pointMap();

    return vtkmesh;
}


// ************************************************************************* //

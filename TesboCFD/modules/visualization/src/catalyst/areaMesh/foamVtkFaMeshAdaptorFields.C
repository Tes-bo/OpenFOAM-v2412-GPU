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

#include "foamVtkFaMeshAdaptor.H"
#include "faMesh.H"
#include "areaFields.H"

// VTK includes
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// Templates (only needed here)
#include "foamVtkFaMeshAdaptorFieldTemplates.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordHashSet Foam::vtk::faMeshAdaptor::knownFields
(
    const wordRes& selectFields
) const
{
    // Quick exit if no volume fields can be converted.
    // This could be refined
    HashTable<wordHashSet> objects = mesh_.thisDb().classes(selectFields);

    wordHashSet allFields(2*objects.size());

    forAllConstIters(objects, iter)
    {
        const word& clsName = iter.key();

        if (fieldTypes::area.contains(clsName))
        {
            allFields += iter.val();
        }
    }

    return allFields;
}


void Foam::vtk::faMeshAdaptor::convertAreaFields
(
    const wordRes& selectFields
)
{
    // Quick exit if no fields can be converted.
    // This could be refined
    HashTable<wordHashSet> objects = mesh_.thisDb().classes(selectFields);

    if (objects.empty())
    {
        return;
    }

    convertAreaFields<scalar>(selectFields);
    convertAreaFields<vector>(selectFields);
    convertAreaFields<sphericalTensor>(selectFields);
    convertAreaFields<symmTensor>(selectFields);
    convertAreaFields<tensor>(selectFields);
}


// ************************************************************************* //

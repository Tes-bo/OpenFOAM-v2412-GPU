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

\*---------------------------------------------------------------------------*/

#include "foamVtkFvMeshAdaptor.H"

// VTK includes
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// Templates (only needed here)
#include "foamVtkFvMeshAdaptorFieldTemplates.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordHashSet Foam::vtk::fvMeshAdaptor::knownFields
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

        if
        (
            fieldTypes::volume.contains(clsName)
            // TODO
            // || fieldTypes::internal.contains(clsName)
        )
        {
            allFields += iter.val();
        }
    }

    return allFields;
}


void Foam::vtk::fvMeshAdaptor::convertVolFields
(
    const wordRes& selectFields
)
{
    // Quick exit if no volume fields can be converted.
    // This could be refined
    HashTable<wordHashSet> objects = mesh_.thisDb().classes(selectFields);

    if (objects.empty())
    {
        return;
    }

    PtrList<patchInterpolator> interpLst;

    if (interpFields_ && patchIds_.size())
    {
        // NOTE: this would be broken with processor patches,
        // but we don't allow them for the catalyst adaptor anyhow

        // patchIds_ are sorted, so the last one is also the max Id

        interpLst.resize(patchIds_.back() + 1);

        for (const label patchId : patchIds_)
        {
            interpLst.set
            (
                patchId,
                new PrimitivePatchInterpolation<primitivePatch>
                (
                    mesh_.boundaryMesh()[patchId]
                )
            );
        }
    }

    convertVolFields<scalar>(interpLst, selectFields);
    convertVolFields<vector>(interpLst, selectFields);
    convertVolFields<sphericalTensor>(interpLst, selectFields);
    convertVolFields<symmTensor>(interpLst, selectFields);
    convertVolFields<tensor>(interpLst, selectFields);

    // TODO
    // convertDimFields<scalar>(interpLst, selectFields);
    // convertDimFields<vector>(interpLst, selectFields);
    // convertDimFields<sphericalTensor>(interpLst, selectFields);
    // convertDimFields<symmTensor>(interpLst, selectFields);
    // convertDimFields<tensor>(interpLst, selectFields);
}


// ************************************************************************* //

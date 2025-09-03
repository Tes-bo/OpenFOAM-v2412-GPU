/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "areaFaMesh.H"
#include "Cloud.H"
#include "IOobjectList.H"
#include "vtkPVFoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// Templates (only needed here)
#include "vtkPVFoamFieldTemplates.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkPVFoam::convertVolFields()
{
    const fvMesh& mesh = *volMeshPtr_;

    const bool interpFields = reader_->GetInterpolateVolFields();
    wordHashSet selectedFields
    (
        getSelected
        (
            reader_->GetVolFieldSelection()
        )
    );

    if (selectedFields.empty())
    {
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects(mesh, dbPtr_().timeName());
    objects.filterKeys(selectedFields);

    if (objects.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        forAllConstIters(objects, iter)
        {
            Info<< "  " << (*iter)->name()
                << " == " << (*iter)->objectPath() << nl;
        }
        printMemory();
    }


    PtrList<patchInterpolator> interpLst;

    if (interpFields)
    {
        interpLst.resize(mesh.boundaryMesh().size());

        forAll(interpLst, i)
        {
            interpLst.set
            (
                i,
                new PrimitivePatchInterpolation<primitivePatch>
                (
                    mesh.boundaryMesh()[i]
                )
            );
        }
    }

    convertVolFields<scalar>(mesh, interpLst, objects);
    convertVolFields<vector>(mesh, interpLst, objects);
    convertVolFields<sphericalTensor>(mesh, interpLst, objects);
    convertVolFields<symmTensor>(mesh, interpLst, objects);
    convertVolFields<tensor>(mesh, interpLst, objects);

    convertDimFields<scalar>(mesh, interpLst, objects);
    convertDimFields<vector>(mesh, interpLst, objects);
    convertDimFields<sphericalTensor>(mesh, interpLst, objects);
    convertDimFields<symmTensor>(mesh, interpLst, objects);
    convertDimFields<tensor>(mesh, interpLst, objects);

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertPointFields()
{
    const fvMesh& mesh = *volMeshPtr_;

    wordHashSet selectedFields
    (
        getSelected
        (
            reader_->GetPointFieldSelection()
        )
    );

    if (selectedFields.empty())
    {
        DebugInfo << "no point fields selected" << nl;
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects(mesh, dbPtr_().timeName());
    objects.filterKeys(selectedFields);

    if (objects.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> convert volume -> point fields" << nl;
        forAllConstIters(objects, iter)
        {
            Info<< "  " << (*iter)->name()
                << " == " << (*iter)->objectPath() << nl;
        }
        printMemory();
    }

    // Construct interpolation on the raw mesh
    const pointMesh& pMesh = pointMesh::New(mesh);

    convertPointFields<scalar>(pMesh, objects);
    convertPointFields<vector>(pMesh, objects);
    convertPointFields<sphericalTensor>(pMesh, objects);
    convertPointFields<symmTensor>(pMesh, objects);
    convertPointFields<tensor>(pMesh, objects);

    if (debug)
    {
        Info<< "<end> convert volume -> point fields" << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertAreaFields()
{
    if (!areaMeshPtr_)
    {
        return;
    }

    const faMesh& mesh = *areaMeshPtr_;

    vtkDataArraySelection* select = reader_->GetVolFieldSelection();

    wordHashSet selectedFields
    (
        getSelected(select)
    );

    if (selectedFields.empty())
    {
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects(mesh.mesh(), dbPtr_().timeName());
    objects.filterKeys(selectedFields);

    if (objects.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        forAllConstIters(objects, iter)
        {
            Info<< "  " << (*iter)->name()
                << " == " << (*iter)->objectPath() << nl;
        }
        printMemory();
    }

    convertAreaFields<scalar>(mesh, objects);
    convertAreaFields<vector>(mesh, objects);
    convertAreaFields<sphericalTensor>(mesh, objects);
    convertAreaFields<symmTensor>(mesh, objects);
    convertAreaFields<tensor>(mesh, objects);

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


void Foam::vtkPVFoam::convertLagrangianFields()
{
    const List<label> partIds =
        rangeClouds_.intersection(selectedPartIds_);

    const fvMesh& mesh = *volMeshPtr_;

    wordHashSet selectedFields
    (
        getSelected
        (
            reader_->GetLagrangianFieldSelection()
        )
    );

    if (selectedFields.empty())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> " << FUNCTION_NAME << nl;
        printMemory();
    }

    for (const auto partId : partIds)
    {
        const auto& longName = selectedPartIds_[partId];
        const word cloudName = getFoamName(longName);

        auto iter = cachedVtp_.find(longName);
        if (!iter.good() || !iter.val().dataset)
        {
            // Should not happen, but for safety require a vtk geometry
            continue;
        }
        auto dataset = iter.val().dataset;

        // Get the Lagrangian fields for this time and this cloud
        // but only keep selected fields
        // the region name is already in the mesh db
        IOobjectList objects
        (
            mesh,
            dbPtr_().timeName(),
            cloud::prefix/cloudName
        );
        objects.filterKeys(selectedFields);

        if (objects.empty())
        {
            continue;
        }

        if (debug)
        {
            Info<< "converting OpenFOAM lagrangian fields" << nl;
            forAllConstIters(objects, iter)
            {
                Info<< "  " << (*iter)->name()
                    << " == " << (*iter)->objectPath() << nl;
            }
        }

        convertLagrangianFields<label>(objects, dataset);
        convertLagrangianFields<scalar>(objects, dataset);
        convertLagrangianFields<vector>(objects, dataset);
        convertLagrangianFields<sphericalTensor>(objects, dataset);
        convertLagrangianFields<symmTensor>(objects, dataset);
        convertLagrangianFields<tensor>(objects, dataset);
    }

    if (debug)
    {
        Info<< "<end> " << FUNCTION_NAME << nl;
        printMemory();
    }
}


// ************************************************************************* //

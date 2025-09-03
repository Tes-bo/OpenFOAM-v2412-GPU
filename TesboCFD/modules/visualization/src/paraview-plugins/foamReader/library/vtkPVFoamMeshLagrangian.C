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
#include "Cloud.H"
#include "fvMesh.H"
#include "IOobjectList.H"
#include "passiveParticle.H"

// VTK includes
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vtkSmartPointer<vtkPolyData> Foam::vtkPVFoam::lagrangianVTKMesh
(
    const polyMesh& mesh,
    const word& cloudName
) const
{
    vtkSmartPointer<vtkPolyData> vtkmesh;

    if (debug)
    {
        Info<< "<beg> lagrangianVTKMesh - timePath "
            << mesh.time().timePath()/cloud::prefix/cloudName << nl;
        printMemory();
    }


    // The region name is already in the mesh db
    IOobjectList cloudObjects
    (
        mesh,
        mesh.time().timeName(),
        cloud::prefix/cloudName
    );

    if
    (
        cloudObjects.contains("positions")
     || cloudObjects.contains("coordinates")
    )
    {
        Cloud<passiveParticle> parcels(mesh, cloudName, false);

        DebugInfo << "cloud with " << parcels.size() << " parcels" << nl;

        auto vtkpoints = vtkSmartPointer<vtkPoints>::New();
        vtkpoints->SetNumberOfPoints(parcels.size());

        vtkIdType particleId = 0;
        for (const passiveParticle& p : parcels)
        {
            vtkpoints->SetPoint(particleId, p.position().v_);
            ++particleId;
        }

        vtkmesh = vtkSmartPointer<vtkPolyData>::New();
        vtkmesh->SetPoints(vtkpoints);
        vtkmesh->SetVerts(vtk::Tools::identityVertices(parcels.size()));
    }

    if (debug)
    {
        Info<< "<end> lagrangianVTKMesh" << nl;
        printMemory();
    }

    return vtkmesh;
}


// ************************************************************************* //

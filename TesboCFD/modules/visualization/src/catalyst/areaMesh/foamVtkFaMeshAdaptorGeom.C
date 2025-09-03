/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::faMeshAdaptor::convertGeometryInternal()
{
    const auto& longName = internalName;

    foamVtpData& vtpData = cachedVtp_(longName);

    vtkSmartPointer<vtkPolyData> vtkgeom;
    if (vtpData.nPoints())
    {
        if (meshState_ == polyMesh::UNCHANGED)
        {
            // Without movement is easy.
            DebugInfo
                << "reuse " << longName << nl;

            vtpData.reuse();
            return;
        }
        else if (meshState_ == polyMesh::POINTS_MOVED)
        {
            DebugInfo
                << "move points " << longName << nl;

            vtkgeom = vtpData.getCopy();
            vtkgeom->SetPoints
            (
                vtk::Tools::Patch::points(mesh_.patch())
            );
        }
    }

    if (!vtkgeom)
    {
        DebugInfo
            << "Nothing usable from cache - create new geometry" << nl;

        vtkgeom = vtk::Tools::Patch::mesh(mesh_.patch());
    }

    vtpData.set(vtkgeom);
}


// ************************************************************************* //

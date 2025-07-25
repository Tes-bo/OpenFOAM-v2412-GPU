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
#include "Pstream.H"

// OpenFOAM includes
#include "faMesh.H"
#include "Time.H"

// VTK includes
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{
    defineTypeNameAndDebug(faMeshAdaptor, 0);
}
} // End namespace Foam

const Foam::word Foam::vtk::faMeshAdaptor::internalName("internal");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::faMeshAdaptor::faMeshAdaptor
(
    const faMesh& mesh
)
:
    mesh_(mesh),
    meshState_(polyMesh::TOPO_CHANGE)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::faMeshAdaptor::updateContent(const wordRes& selectFields)
{
    // Update cached, saved, unneed values.
    if
    (
        meshState_ == polyMesh::TOPO_CHANGE
     || meshState_ == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        // Eliminate cached values that would be unreliable
        forAllIters(cachedVtp_, iter)
        {
            iter.val().clearGeom();
            iter.val().clear();
        }
    }

    convertGeometryInternal();
    convertAreaFields(selectFields);
    meshState_ = polyMesh::UNCHANGED;
}


vtkSmartPointer<vtkMultiPieceDataSet>
Foam::vtk::faMeshAdaptor::output(const wordRes& select)
{
    updateContent(select);

    // All individual datasets are vtkMultiPieceDataSet for improved
    // handling downstream.

    label rank = 0;
    label nproc = 1;

    if (UPstream::parRun())
    {
        rank  = UPstream::myProcNo();
        nproc = UPstream::nProcs();
    }

    auto output = vtkSmartPointer<vtkMultiPieceDataSet>::New();

    vtkSmartPointer<vtkPolyData> vtkmesh;

    // MESH
    {
        const auto& longName = internalName;
        auto iter = cachedVtp_.find(longName);
        if (iter.good() && iter.val().dataset)
        {
            vtkmesh = iter.val().dataset;
        }
    }

    output->SetNumberOfPieces(nproc);
    output->SetPiece(rank, vtkmesh);

    return output;
}


void Foam::vtk::faMeshAdaptor::updateState(polyMesh::readUpdateState state)
{
    // Only move to worse states
    switch (state)
    {
        case polyMesh::UNCHANGED:
            break;

        case polyMesh::POINTS_MOVED:
            if (meshState_ == polyMesh::UNCHANGED)
            {
                meshState_ = polyMesh::POINTS_MOVED;
            }
            break;

        case polyMesh::TOPO_CHANGE:
        case polyMesh::TOPO_PATCH_CHANGE:
            meshState_ = polyMesh::TOPO_CHANGE;
            break;
    }
}


// ************************************************************************* //

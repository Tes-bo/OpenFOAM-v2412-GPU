/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

// OpenFOAM includes
#include "surface.H"
#include "runTimePostProcessing.H"

#include "foamVtkTools.H"
#include "polySurface.H"

// VTK includes
#include "vtkMultiPieceDataSet.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

vtkSmartPointer<vtkMultiPieceDataSet>
Foam::functionObjects::runTimePostPro::surface::gatherSurfacePieces
(
    const polySurface* surf
) const
{
    auto multiPiece = vtkSmartPointer<vtkMultiPieceDataSet>::New();
    multiPiece->SetNumberOfPieces(Pstream::nProcs());

    if (!needsCollective())
    {
        // Simple case (serial-serial, parallel-parallel)

        if (surf)
        {
            multiPiece->SetPiece
            (
                Pstream::myProcNo(),
                Foam::vtk::Tools::Patch::mesh(*surf)
            );
        }
    }
    else if (Pstream::master())
    {
        // Gather pieces on master

        if (surf)
        {
            // Add myself

            multiPiece->SetPiece
            (
                Pstream::myProcNo(),
                Foam::vtk::Tools::Patch::mesh(*surf)
            );
        }

        // Receive surfaces
        pointField points;
        faceList faces;

        for (const int proci : Pstream::subProcs())
        {
            IPstream fromProc(Pstream::commsTypes::scheduled, proci);

            points.clear();
            faces.clear();

            fromProc >> points >> faces;

            if (points.size())
            {
                multiPiece->SetPiece
                (
                    proci,
                    Foam::vtk::Tools::Patch::mesh(points, faces)
                );
            }
        }
    }
    else
    {
        // Sub-process - send surfaces

        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        if (surf)
        {
            toMaster << surf->points() << surf->faces();
        }
        else
        {
            toMaster << pointField() << faceList();
        }
    }

    return multiPiece;
}


vtkSmartPointer<vtkMultiPieceDataSet>
Foam::functionObjects::runTimePostPro::surface::gatherFaceCentres
(
    const polySurface* surf
) const
{
    auto multiPiece = vtkSmartPointer<vtkMultiPieceDataSet>::New();
    multiPiece->SetNumberOfPieces(Pstream::nProcs());

    if (!needsCollective())
    {
        // Simple case

        if (surf)
        {
            auto dataset = vtkSmartPointer<vtkPolyData>::New();

            auto geom = vtkSmartPointer<vtkPolyData>::New();

            geom->SetPoints(Foam::vtk::Tools::Patch::faceCentres(*surf));
            geom->SetVerts(Foam::vtk::Tools::identityVertices(surf->nFaces()));

            multiPiece->SetPiece(Pstream::myProcNo(), geom);
        }
    }
    else if (Pstream::master())
    {
        // Gather pieces (face centres) on master

        if (surf)
        {
            // Add myself

            auto geom = vtkSmartPointer<vtkPolyData>::New();

            geom->SetPoints(Foam::vtk::Tools::Patch::faceCentres(*surf));
            geom->SetVerts(Foam::vtk::Tools::identityVertices(surf->nFaces()));

            multiPiece->SetPiece(Pstream::myProcNo(), geom);
        }

        // Receive points
        pointField points;

        for (const int proci : Pstream::subProcs())
        {
            IPstream fromProc(Pstream::commsTypes::scheduled, proci);

            points.clear();

            fromProc >> points;

            if (points.size())
            {
                multiPiece->SetPiece
                (
                    proci,
                    Foam::vtk::Tools::Vertices(points)
                );
            }
        }
    }
    else
    {
        // Sub-process - send face centres

        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        if (surf)
        {
            toMaster << surf->faceCentres();
        }
        else
        {
            toMaster << pointField();
        }
    }

    return multiPiece;
}


// ************************************************************************* //

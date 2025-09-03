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
#include "geometryCloud.H"
#include "cloud.H"
#include "runTimePostProcessing.H"

#include "foamVtkTools.H"

// VTK includes
#include "vtkMultiPieceDataSet.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

vtkSmartPointer<vtkMultiPieceDataSet>
Foam::functionObjects::runTimePostPro::geometryCloud::gatherCloud
(
    const objectRegistry& obrTmp
) const
{
    auto multiPiece = vtkSmartPointer<vtkMultiPieceDataSet>::New();
    multiPiece->SetNumberOfPieces(Pstream::nProcs());

    const auto* pointsPtr = cloud::findIOPosition(obrTmp);

    if (!needsCollective())
    {
        // Simple case (serial-serial, parallel-parallel)

        if (pointsPtr && pointsPtr->size())
        {
            multiPiece->SetPiece
            (
                Pstream::myProcNo(),
                Foam::vtk::Tools::Vertices(*pointsPtr)
            );
        }
    }
    else if (Pstream::master())
    {
        // Gather pieces on master

        if (pointsPtr && pointsPtr->size())
        {
            // Add myself

            multiPiece->SetPiece
            (
                Pstream::myProcNo(),
                Foam::vtk::Tools::Vertices(*pointsPtr)
            );
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
        // Sub-process - send points

        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        if (pointsPtr)
        {
            toMaster << *pointsPtr;
        }
        else
        {
            toMaster << pointField();
        }
    }

    return multiPiece;
}


// ************************************************************************* //

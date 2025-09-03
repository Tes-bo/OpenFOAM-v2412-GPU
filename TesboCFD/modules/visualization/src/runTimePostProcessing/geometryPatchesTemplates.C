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
#include "geometryPatches.H"
#include "fvMesh.H"
#include "volFields.H"
#include "foamVtkTools.H"

// VTK includes
#include "vtkCellData.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkPointData.h"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
int Foam::functionObjects::runTimePostPro::geometryPatches::addPatchField
(
    vtkMultiPieceDataSet* multiPiece,
    const labelListList& patchIds,
    const GeometricField<Type, fvPatchField, volMesh>* fldptr,
    const word& fieldName
) const
{
    if (!multiPiece || !fldptr)
    {
        return 0;
    }

    const int nCmpt(pTraits<Type>::nComponents);

    const auto& bf = fldptr->boundaryField();

    label pieceId = 0;

    if (!needsCollective())
    {
        // Simple case (serial-serial, parallel-parallel)

        for (int proci=0; proci < Pstream::myProcNo(); ++proci)
        {
            pieceId += patchIds[proci].size();
        }

        for (const label patchId : patchIds[Pstream::myProcNo()])
        {
            const auto& pf = bf[patchId];

            auto vtkfield =
            (
                nearCellValue_
              ? Foam::vtk::Tools::convertFieldToVTK<Type>
                (
                    fieldName,
                    pf.patchInternalField()()
                )
              : Foam::vtk::Tools::convertFieldToVTK<Type>
                (
                    fieldName,
                    pf
                )
            );

            auto piece = multiPiece->GetPiece(pieceId);

            if (piece->GetNumberOfCells() == piece->GetNumberOfPoints())
            {
                // Only has verts
                piece->GetPointData()->AddArray(vtkfield);
            }
            else
            {
                piece->GetCellData()->AddArray(vtkfield);
            }

            ++pieceId;
        }
    }
    else if (Pstream::master())
    {
        // Gather pieces on master

        // Add myself
        for (const label patchId : patchIds[Pstream::myProcNo()])
        {
            const auto& pf = bf[patchId];

            auto vtkfield =
            (
                nearCellValue_
              ? Foam::vtk::Tools::convertFieldToVTK<Type>
                (
                    fieldName,
                    pf.patchInternalField()()
                )
              : Foam::vtk::Tools::convertFieldToVTK<Type>
                (
                    fieldName,
                    pf
                )
            );

            auto piece = multiPiece->GetPiece(pieceId);

            if (piece->GetNumberOfCells() == piece->GetNumberOfPoints())
            {
                // Only has verts
                piece->GetPointData()->AddArray(vtkfield);
            }
            else
            {
                piece->GetCellData()->AddArray(vtkfield);
            }

            ++pieceId;
        }

        // Receive field
        Field<Type> recv;

        for (const int proci : Pstream::subProcs())
        {
            const label nPatches = patchIds[proci].size();

            if (!nPatches)
            {
                continue;
            }

            IPstream fromProc(Pstream::commsTypes::scheduled, proci);

            for (label recvi=0; recvi < nPatches; ++recvi)
            {
                recv.clear();

                fromProc >> recv;

                auto vtkfield = Foam::vtk::Tools::convertFieldToVTK<Type>
                (
                    fieldName,
                    recv
                );

                auto piece = multiPiece->GetPiece(pieceId);

                if (piece->GetNumberOfCells() == piece->GetNumberOfPoints())
                {
                    // Only has verts
                    piece->GetPointData()->AddArray(vtkfield);
                }
                else
                {
                    piece->GetCellData()->AddArray(vtkfield);
                }

                ++pieceId;
            }
        }
    }
    else
    {
        // Sub-process - send fields

        if (patchIds[Pstream::myProcNo()].size())
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            for (const label patchId : patchIds[Pstream::myProcNo()])
            {
                const auto& pf = bf[patchId];

                if (nearCellValue_)
                {
                    // Send patch field internal values
                    toMaster << pf.patchInternalField()();
                }
                else
                {
                    // Send patch field values
                    toMaster << static_cast<const Field<Type>&>(pf);
                }
            }
        }
    }

    return nCmpt;
}


// ************************************************************************* //

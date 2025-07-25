/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
#include "volumeFilter.H"
#include "fvMesh.H"
#include "volFields.H"
#include "foamVtkTools.H"

// VTK includes
#include "vtkCellData.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkPointData.h"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::runTimePostPro::volumeFilter::addDimField
(
    vtkDataSet* piece,
    const vtk::vtuAdaptor& adaptor,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    vtkSmartPointer<vtkFloatArray> vtkdata;

    const auto* dimptr =
        dynamic_cast<const DimensionedField<Type, volMesh>*>(ioptr);

    if (dimptr && !vtkdata)
    {
        vtkdata = adaptor.convertField(*dimptr);
    }

    const auto* volptr =
        dynamic_cast<const GeometricField<Type, fvPatchField, volMesh>*>(ioptr);

    if (volptr && !vtkdata)
    {
        vtkdata = adaptor.convertField(volptr->internalField());
    }

    if (vtkdata)
    {
        piece->GetCellData()->AddArray(vtkdata);
        return true;
    }

    return false;
}


template<class Type>
int Foam::functionObjects::runTimePostPro::volumeFilter::addDimField
(
    vtkMultiPieceDataSet* multiPiece,
    const vtk::vtuAdaptor& adaptor,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    if (!multiPiece)
    {
        return 0;
    }

    const int nCmpt(pTraits<Type>::nComponents);

    if (!needsCollective())
    {
        // Simple case (serial-serial, parallel-parallel)

        auto piece = multiPiece->GetPiece(Pstream::myProcNo());

        if
        (
            addDimField<Type>
            (
                piece,
                adaptor,
                ioptr,
                fieldName
            )
        )
        {
            return nCmpt;
        }
    }
    else
    {
    }

    return 0;
}


// ************************************************************************* //

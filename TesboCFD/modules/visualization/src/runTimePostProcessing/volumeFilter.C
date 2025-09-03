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
#include "runTimePostProcessing.H"

// VTK includes
#include "vtkCellData.h"
#include "vtkCellDataToPointData.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkCompositeDataSet.h"
#include "vtkCompositePolyDataMapper.h"
#include "vtkMultiPieceDataSet.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::runTimePostPro::volumeFilter::volumeFilter
(
    const runTimePostProcessing& parent,
    const dictionary& dict,
    const HashPtrTable<Function1<vector>>& colours
)
:
    surface(parent, dict, colours)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

vtkSmartPointer<vtkMultiPieceDataSet>
Foam::functionObjects::runTimePostPro::volumeFilter::mesh
(
    Foam::vtk::vtuAdaptor& adaptor
) const
{
    auto multiPiece = vtkSmartPointer<vtkMultiPieceDataSet>::New();

    multiPiece->SetNumberOfPieces(Pstream::nProcs());
    multiPiece->SetPiece
    (
        Pstream::myProcNo(),
        adaptor.internal(parent().mesh())
    );

    return multiPiece;
}


bool Foam::functionObjects::runTimePostPro::volumeFilter::addDimField
(
    vtkDataSet* piece,
    const vtk::vtuAdaptor& adaptor,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    return (piece && ioptr) &&
    (
        addDimField<scalar>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<vector>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<sphericalTensor>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<symmTensor>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<tensor>
        (
            piece, adaptor, ioptr, fieldName
        )
    );
}


int Foam::functionObjects::runTimePostPro::volumeFilter::addDimField
(
    vtkMultiPieceDataSet* piece,
    const vtk::vtuAdaptor& adaptor,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    return (piece && ioptr) &&
    (
        addDimField<scalar>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<vector>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<sphericalTensor>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<symmTensor>
        (
            piece, adaptor, ioptr, fieldName
        )
     || addDimField<tensor>
        (
            piece, adaptor, ioptr, fieldName
        )
    );
}


// ************************************************************************* //

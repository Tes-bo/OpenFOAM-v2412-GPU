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

#include "foamVtkTools.H"
#include "polySurfaceFields.H"
#include "polySurfacePointFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMeshType>
bool Foam::functionObjects::runTimePostPro::surface::addField
(
    vtkDataSet* piece,
    const Field<Type>& fld,
    const word& fieldName
) const
{
    if (!piece) return false;

    auto vtkfield = Foam::vtk::Tools::convertFieldToVTK<Type>(fieldName, fld);

    if (piece->GetNumberOfCells() == piece->GetNumberOfPoints())
    {
        // Only has verts
        piece->GetPointData()->AddArray(vtkfield);
    }
    else
    {
        Foam::vtk::Tools::FieldAccess<GeoMeshType>()(piece)->AddArray(vtkfield);
    }

    return true;
}


template<class Type, class GeoMeshType>
bool Foam::functionObjects::runTimePostPro::surface::addDimField
(
    vtkDataSet* piece,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    const auto* fldptr =
        dynamic_cast<const DimensionedField<Type, GeoMeshType>*>(ioptr);

    if (fldptr)
    {
        return addField<Type, GeoMeshType>
        (
            piece,
            fldptr->field(),
            fieldName
        );
    }

    return false;
}


template<class GeoMeshType>
bool Foam::functionObjects::runTimePostPro::surface::addDimField
(
    vtkDataSet* piece,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    return (piece && ioptr) &&
    (
        addDimField<scalar, GeoMeshType>
        (
            piece, ioptr, fieldName
        )
     || addDimField<vector, GeoMeshType>
        (
            piece, ioptr, fieldName
        )
     || addDimField<sphericalTensor, GeoMeshType>
        (
            piece, ioptr, fieldName
        )
     || addDimField<symmTensor, GeoMeshType>
        (
            piece, ioptr, fieldName
        )
     || addDimField<tensor, GeoMeshType>
        (
            piece, ioptr, fieldName
        )
    );
}


template<class Type, class GeoMeshType>
bool Foam::functionObjects::runTimePostPro::surface::addDimField
(
    vtkMultiPieceDataSet* multiPiece,
    const DimensionedField<Type, GeoMeshType>* fldptr,
    const word& fieldName
) const
{
    if (!multiPiece)
    {
        return false;
    }

    if (!needsCollective())
    {
        // Simple case (serial-serial, parallel-parallel)

        return fldptr &&
            addField<Type, GeoMeshType>
            (
                multiPiece->GetPiece(Pstream::myProcNo()),
                fldptr->field(),
                fieldName
            );
    }


    // Gather fields
    const bool ok = returnReduce((fldptr != nullptr), orOp<bool>());

    if (!ok)
    {
        return false;
    }

    if (Pstream::master())
    {
        if (fldptr)
        {
            // My field data
            addField<Type, GeoMeshType>
            (
                multiPiece->GetPiece(Pstream::myProcNo()),
                fldptr->field(),
                fieldName
            );
        }

        // Receive field data
        Field<Type> recv;

        for (const int proci : Pstream::subProcs())
        {
            IPstream fromProc(Pstream::commsTypes::scheduled, proci);

            recv.clear();

            fromProc >> recv;

            if (recv.size())
            {
                addField<Type, GeoMeshType>
                (
                    multiPiece->GetPiece(proci),
                    recv,
                    fieldName
                );
            }
        }
    }
    else
    {
        // Sub-process - send field data

        OPstream toMaster
        (
            Pstream::commsTypes::scheduled,
            Pstream::masterNo()
        );

        if (fldptr)
        {
            toMaster << fldptr->field();
        }
        else
        {
            toMaster << Field<Type>();
        }
    }

    return ok;
}


template<class Type, class GeoMeshType>
bool Foam::functionObjects::runTimePostPro::surface::addDimField
(
    vtkMultiPieceDataSet* multiPiece,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    return addDimField<Type, GeoMeshType>
    (
        multiPiece,
        dynamic_cast<const DimensionedField<Type, GeoMeshType>*>(ioptr),
        fieldName
    );
}


template<class GeoMeshType>
bool Foam::functionObjects::runTimePostPro::surface::addDimField
(
    vtkMultiPieceDataSet* multiPiece,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    return (multiPiece) &&
    (
        addDimField<scalar, GeoMeshType>
        (
            multiPiece, ioptr, fieldName
        )
     || addDimField<vector, GeoMeshType>
        (
            multiPiece, ioptr, fieldName
        )
     || addDimField<sphericalTensor, GeoMeshType>
        (
            multiPiece, ioptr, fieldName
        )
     || addDimField<symmTensor, GeoMeshType>
        (
            multiPiece, ioptr, fieldName
        )
     || addDimField<tensor, GeoMeshType>
        (
            multiPiece, ioptr, fieldName
        )
    );
}


template<class GeoMeshType>
bool Foam::functionObjects::runTimePostPro::surface::addDimField
(
    vtkMultiPieceDataSet* multiPiece,
    const polySurface* surf,
    const word& fieldName
) const
{
    const regIOobject* ioptr =
    (
        surf
      ? surf->findFieldObject<GeoMeshType>(fieldName)
      : nullptr
    );

    return addDimField<GeoMeshType>(multiPiece, ioptr, fieldName);
}


// ************************************************************************* //

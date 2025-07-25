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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::runTimePostPro::geometryCloud::addField
(
    vtkDataSet* piece,
    const Field<Type>& fld,
    const word& fieldName
) const
{
    if (!piece) return false;

    auto vtkfield = Foam::vtk::Tools::convertFieldToVTK<Type>(fieldName, fld);

    // Only has verts
    piece->GetPointData()->AddArray(vtkfield);

    return true;
}


template<class Type>
bool Foam::functionObjects::runTimePostPro::geometryCloud::addCloudField
(
    vtkMultiPieceDataSet* multiPiece,
    const IOField<Type>* fldptr,
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
            addField<Type>
            (
                multiPiece->GetPiece(Pstream::myProcNo()),
                *fldptr,
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
            addField<Type>
            (
                multiPiece->GetPiece(Pstream::myProcNo()),
                *fldptr,
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
                addField<Type>
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
            toMaster << *fldptr;
        }
        else
        {
            toMaster << Field<Type>();
        }
    }

    return ok;
}


template<class Type>
bool Foam::functionObjects::runTimePostPro::geometryCloud::addCloudField
(
    vtkMultiPieceDataSet* multiPiece,
    const regIOobject* ioptr,
    const word& fieldName
) const
{
    return addCloudField<Type>
    (
        multiPiece,
        dynamic_cast<const IOField<Type>*>(ioptr),
        fieldName
    );
}


// ************************************************************************* //

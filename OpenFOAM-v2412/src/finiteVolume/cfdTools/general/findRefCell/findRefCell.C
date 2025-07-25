/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "findRefCell.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::setRefCell
(
    const volScalarField& field,
    const volScalarField& fieldRef,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    if (fieldRef.needReference() || forceReference)
    {
        const word refCellName = field.name() + "RefCell";
        const word refPointName = field.name() + "RefPoint";
        const word refValueName = field.name() + "RefValue";

        if (dict.found(refCellName))
        {
            if (Pstream::master())
            {
                dict.readEntry(refCellName, refCelli);

                if (refCelli < 0 || refCelli >= field.mesh().nCells())
                {
                    FatalIOErrorInFunction(dict)
                        << "Illegal master cellID " << refCelli
                        << ". Should be 0.." << field.mesh().nCells()
                        << exit(FatalIOError);
                }
            }
            else
            {
                refCelli = -1;
            }
        }
        else if (dict.found(refPointName))
        {
            point refPointi(dict.get<point>(refPointName));

            // Try fast approximate search avoiding octree construction
            refCelli = field.mesh().findCell(refPointi, polyMesh::FACE_PLANES);

            label hasRef = (refCelli >= 0 ? 1 : 0);
            label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

            // If reference cell not found, use octree search
            // with cell tet-decompositoin
            if (sumHasRef != 1)
            {
                refCelli = field.mesh().findCell(refPointi);

                hasRef = (refCelli >= 0 ? 1 : 0);
                sumHasRef = returnReduce<label>(hasRef, sumOp<label>());
            }

            if (sumHasRef != 1)
            {
                FatalIOErrorInFunction(dict)
                    << "Unable to set reference cell for field " << field.name()
                    << nl << "    Reference point " << refPointName
                    << " " << refPointi
                    << " found on " << sumHasRef << " domains (should be one)"
                    << nl << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Unable to set reference cell for field " << field.name()
                << nl
                << "    Please supply either " << refCellName
                << " or " << refPointName << nl << exit(FatalIOError);
        }

        dict.readEntry(refValueName, refValue);

        return true;
    }
    else
    {
        refCelli = -1;
    }

    return false;
}


bool Foam::setRefCell
(
    const volScalarField& field,
    const dictionary& dict,
    label& refCelli,
    scalar& refValue,
    const bool forceReference
)
{
    return setRefCell(field, field, dict, refCelli, refValue, forceReference);
}


Foam::scalar Foam::getRefCellValue
(
    const volScalarField& field,
    const label refCelli
)
{
    #ifdef FULLDEBUG
    if (refCelli >= field.mesh().nCells())
    {
        FatalErrorInFunction
            << "Illegal reference cellID " << refCelli
            << ". Mesh has " << field.mesh().nCells() << ". cells."
            << exit(FatalError);
    }
    #endif
    scalar refCellValue = (refCelli >= 0 ? field[refCelli] : 0.0);

    // Currently distributing the value to all processors. This is generally
    // not needed since only the processor holding the reference cell needs
    // it. Tdb.
    return returnReduce(refCellValue, sumOp<scalar>());
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "faMeshBoundaryHalo.H"
#include "IndirectList.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::faMeshBoundaryHalo::distributeSparse
(
    List<Type>& fld,
    const labelUList& sparseInputLocations,
    const labelUList& compactOutputMapping
) const
{
    if (!Pstream::parRun())
    {
        return;  // No-op in serial
    }

    // Construct data in compact addressing
    List<Type> compactFld(constructSize(), Zero);

    if (sparseInputLocations.empty())
    {
        // Copy in as dense field
        forAll(fld, i)
        {
            compactFld[i] = fld[i];
        }
    }
    else
    {
        if (fld.size() != sparseInputLocations.size())
        {
            FatalErrorInFunction
                << "Input field size (" << fld.size()
                << " != sparse ids size ("
                << sparseInputLocations.size() << ")\n"
                << exit(FatalError);
        }

        // Copy in sparse locations
        forAll(sparseInputLocations, i)
        {
            const label idx = sparseInputLocations[i];
            if (idx != -1)
            {
                compactFld[idx] = fld[i];
            }
        }
    }


    // Pull all data
    mapDistributeBase::distribute<Type>(compactFld);

    // Rewrite to output
    fld = UIndirectList<Type>(compactFld, compactOutputMapping);
}


template<class Type>
void Foam::faMeshBoundaryHalo::distributeSparse
(
    List<Type>& fld,
    const labelUList& sparseInputLocations
) const
{
    this->distributeSparse(fld, sparseInputLocations, boundaryToCompact_);
}


template<class Type>
void Foam::faMeshBoundaryHalo::distributeSparse(List<Type>& fld) const
{
    this->distributeSparse(fld, inputMeshFaces_, boundaryToCompact_);
}


// ************************************************************************* //

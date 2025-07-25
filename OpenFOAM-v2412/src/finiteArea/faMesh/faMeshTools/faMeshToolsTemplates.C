/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "edgeFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::faMeshTools::flattenEdgeField
(
    const GeometricField<Type, faePatchField, edgeMesh>& fld,
    const bool primitiveOrdering
)
{
    const faMesh& mesh = fld.mesh();

    auto tresult = tmp<Field<Type>>::New(mesh.nEdges(), Zero);
    auto& result = tresult.ref();

    // Internal field
    result.slice(0, fld.size()) = fld.primitiveField();

    label start = fld.size();

    // Boundary fields
    forAll(fld.boundaryField(), patchi)
    {
        const labelList& edgeLabels = mesh.boundary()[patchi].edgeLabels();
        const label len = edgeLabels.size();
        const auto& pfld = fld.boundaryField()[patchi];

        // Only assign when field size matches underlying patch size
        // ie, skip 'empty' patches etc

        if (len == pfld.size())
        {
            if (primitiveOrdering)
            {
                // In primitive patch order
                UIndirectList<Type>(result, edgeLabels) = pfld;
            }
            else
            {
                // In sub-list (slice) order
                result.slice(start, len) = pfld;
            }
        }

        start += len;
    }

    return tresult;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurfaceTopo::interpolateTemplate
(
    const Field<Type>& cellData,
    const Field<Type>& pointData
) const
{
    auto tfld = tmp<Field<Type>>::New(pointToVerts_.size());
    auto& fld = tfld.ref();

    forAll(pointToVerts_, i)
    {
        const edge& verts = pointToVerts_[i];
        Type& val = fld[i];

        scalar s0;
        Type v0;
        {
            label idx = verts.first();
            if (idx < mesh_.nPoints())
            {
                // Point index
                s0 = pVals_[idx];
                v0 = pointData[idx];
            }
            else
            {
                // Cell index
                idx -= mesh_.nPoints();
                s0 = cVals_[idx];
                v0 = cellData[idx];
            }
        }

        scalar s1;
        Type v1;
        {
            label idx = verts.second();
            if (idx == verts.first())
            {
                // Duplicate index (ie, snapped)
                val = v0;
                continue;
            }
            else if (idx < mesh_.nPoints())
            {
                // Point index
                s1 = pVals_[idx];
                v1 = pointData[idx];
            }
            else
            {
                // Cell index
                idx -= mesh_.nPoints();
                s1 = cVals_[idx];
                v1 = cellData[idx];
            }
        }

        const scalar d = s1-s0;
        if (mag(d) > VSMALL)
        {
            const scalar s = (iso_-s0)/d;
            val = s*v1+(1.0-s)*v0;
        }
        else
        {
            val = 0.5*(v0+v1);
        }
    }

    return tfld;
}


// ************************************************************************* //

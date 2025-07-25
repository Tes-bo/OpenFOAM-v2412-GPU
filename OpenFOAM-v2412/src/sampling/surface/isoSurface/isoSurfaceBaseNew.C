/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "isoSurfaceBase.H"
#include "isoSurfaceCell.H"
#include "isoSurfacePoint.H"
#include "isoSurfaceTopo.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::isoSurfaceBase>
Foam::isoSurfaceBase::New
(
    const isoSurfaceParams& params,
    const volScalarField& cellValues,
    const scalarField& pointValues,
    const scalar iso,
    const bitSet& ignoreCells
)
{
    autoPtr<isoSurfaceBase> ptr;

    if (params.algorithm() == isoSurfaceParams::ALGO_POINT)
    {
        ptr.reset
        (
            new isoSurfacePoint
            (
                /* fvMesh implicit from cellValues */
                cellValues,
                pointValues,
                iso,
                params,
                ignoreCells // unused
            )
        );
    }
    else if (params.algorithm() == isoSurfaceParams::ALGO_CELL)
    {
        ptr.reset
        (
            new isoSurfaceCell
            (
                cellValues.mesh(),  // polyMesh
                cellValues.primitiveField(),
                pointValues,
                iso,
                params,
                ignoreCells
            )
        );
    }
    else
    {
        // ALGO_TOPO, ALGO_DEFAULT
        ptr.reset
        (
            new isoSurfaceTopo
            (
                cellValues.mesh(),  // polyMesh
                cellValues.primitiveField(),
                pointValues,
                iso,
                params,
                ignoreCells
            )
        );
    }

    // Cannot really fail (with current logic)
    // if (!ptr)
    // {
    //     FatalErrorInFunction
    //         << "Unknown iso-surface algorithm ("
    //         << int(params.algorithm()) << ")\n"
    //         << exit(FatalError);
    // }

    return ptr;
}


// ************************************************************************* //

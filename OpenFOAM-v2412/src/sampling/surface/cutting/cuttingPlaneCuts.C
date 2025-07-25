/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "cuttingPlane.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Check for face/plane intersection based on crossings
    // Took (-1,0,+1) from plane::sign and packed as (0,1,2).
    // Now use for left shift to obtain (1,2,4).
    //
    // Test accumulated value for an intersection with the plane.
    inline bool intersectsFace
    (
        const PackedList<2>& sides,
        const labelUList& indices
    )
    {
        unsigned accum = 0u;

        for (const label pointi : indices)
        {
            accum |= (1u << sides[pointi]);
        }

        // Accumulated value 3,5,6,7 indicates an intersection
        return (accum == 3 || accum >= 5);
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::cuttingPlane::calcCellCuts
(
    const primitiveMesh& mesh,
    PackedList<2>& sides,
    bitSet& cellCuts
)
{
    const faceList& faces = mesh.faces();
    const pointField& pts = mesh.points();

    const label nCells = mesh.nCells();
    const label nFaces = mesh.nFaces();
    const label nInternFaces = mesh.nInternalFaces();


    // Classify sides of plane (0=BACK, 1=ONPLANE, 2=FRONT) for each point
    {
        const plane& pln = *this;
        const label len = pts.size();

        sides.resize(len);

        // From signed (-1,0,+1) to (0,1,2) for PackedList
        for (label i=0; i < len; ++i)
        {
            sides.set(i, unsigned(1 + pln.sign(pts[i], SMALL)));
        }
    }


    // Detect cells cuts from the face cuts

    label nFaceCuts = 0;

    // 1st face cut of cell
    bitSet hasCut1(nCells);

    // 2nd face cut of cell
    bitSet hasCut2(nCells);

    for (label facei = 0; facei < nInternFaces; ++facei)
    {
        if (intersectsFace(sides, faces[facei]))
        {
            const label own = mesh.faceOwner()[facei];
            const label nei = mesh.faceNeighbour()[facei];

            ++nFaceCuts;

            if (!hasCut1.set(own))
            {
                hasCut2.set(own);
            }
            if (!hasCut1.set(nei))
            {
                hasCut2.set(nei);
            }
        }
    }

    for (label facei = nInternFaces; facei < nFaces; ++facei)
    {
        if (intersectsFace(sides, faces[facei]))
        {
            const label own = mesh.faceOwner()[facei];

            ++nFaceCuts;

            if (!hasCut1.set(own))
            {
                hasCut2.set(own);
            }
        }
    }

    hasCut1.clearStorage();   // Not needed now

    if (cellCuts.size())
    {
        cellCuts.resize(nCells);    // safety
        cellCuts &= hasCut2;        // restrict to cell subset

        if (debug)
        {
            Pout<<"detected " << cellCuts.count() << "/" << nCells
                << " cells cut, subsetted from "
                << hasCut2.count() << "/" << nCells << " cells." << endl;
        }
    }
    else
    {
        cellCuts = std::move(hasCut2);

        if (debug)
        {
            Pout<<"detected " << cellCuts.count() << "/" << nCells
                << " cells cut." << endl;
        }
    }


    const fvMesh* fvMeshPtr = nullptr;
    if (debug && (fvMeshPtr = isA<fvMesh>(mesh)) != nullptr)
    {
        auto tcellCutsDebug = volScalarField::New
        (
            "cuttingPlane.cellCuts",
            IOobjectOption::NO_REGISTER,
            *fvMeshPtr,
            dimensionedScalar(dimless, Foam::zero{})
        );
        auto& cellCutsDebug = tcellCutsDebug.ref();

        auto& fld = cellCutsDebug.primitiveFieldRef();

        for (const label celli : cellCuts)
        {
            fld[celli] = 1;
        }

        Pout<< "Writing cut types:" << cellCutsDebug.objectPath() << endl;
        cellCutsDebug.write();
    }


    return nFaceCuts;
}


// ************************************************************************* //

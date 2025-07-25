/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "primitiveMesh.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::edgeFaces() const
{
    if (!efPtr_)
    {
        if (debug)
        {
            Pout<< "primitiveMesh::edgeFaces() : calculating edgeFaces" << endl;

            if (debug == -1)
            {
                // For checking calls:abort so we can quickly hunt down
                // origin of call
                FatalErrorInFunction
                    << abort(FatalError);
            }
        }

        // Invert faceEdges
        efPtr_ = new labelListList(nEdges());
        invertManyToMany(nEdges(), faceEdges(), *efPtr_);
    }

    return *efPtr_;
}


const Foam::labelList& Foam::primitiveMesh::edgeFaces
(
    const label edgei,
    DynamicList<label>& storage
) const
{
    if (hasEdgeFaces())
    {
        return edgeFaces()[edgei];
    }
    else
    {
        // Use the fact that pointFaces are sorted in incrementing edge order
        // (since they get constructed by inverting the faces which walks
        //  in increasing face order)
        const edge& e = edges()[edgei];
        const labelList& pFaces0 = pointFaces()[e[0]];
        const labelList& pFaces1 = pointFaces()[e[1]];

        label i0 = 0;
        label i1 = 0;

        storage.clear();

        while (i0 < pFaces0.size() && i1 < pFaces1.size())
        {
            const label f0 = pFaces0[i0];
            const label f1 = pFaces1[i1];

            if (f0 < f1)
            {
                ++i0;
            }
            else if (f0 > f1)
            {
                ++i1;
            }
            else
            {
                // Equal face. Check if indeed on consecutive vertices on both
                // faces since it could be that there is an 'ear' where one
                // side is a triangle  and the other side is part of a bigger
                // face (e.g. quad). Now all vertex-vertex pairs on the
                // triangle are edges but there is no cross connection on the
                // bigger face.
                const face& f = faces()[f0];
                const label fp0 = f.find(e[0]);

                if (f[f.fcIndex(fp0)] == e[1] || f[f.rcIndex(fp0)] == e[1])
                {
                    storage.append(f0);
                    ++i0;
                    ++i1;
                }
                else
                {
                    ++i1;
                }
            }
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::edgeFaces(const label edgei) const
{
    return edgeFaces(edgei, labels_);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "blockMesh.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::blockMesh::checkDegenerate() const
{
    const blockList& blocks = *this;

    for (const block& blk : blocks)
    {
        const cellShape& shape = blk.blockShape();

        if (shape.model().index() == cellModel::HEX)
        {
            // Check for collapsed edges
            // - limit to HEX only for now.
            for (label edgei = 0; edgei < shape.nEdges(); ++edgei)
            {
                edge e(shape.edge(edgei));

                if (!e.good())
                {
                    return true;  // Looks like a collapsed edge
                }
            }
        }
    }

    return false;
}


void Foam::blockMesh::check
(
    const polyMesh& bm,
    const dictionary& dict
) const
{
    Info<< nl << "Check topology" << endl;

    bool ok = true;

    // Check for duplicate curved edge definitions
    forAll(edges_, cei)
    {
        for (label cej=cei+1; cej<edges_.size(); cej++)
        {
            if (edges_[cei].compare(edges_[cej]) != 0)
            {
                Info<< "    Curved edge ";
                edges_[cej].write(Info, dict);
                Info<< "    is a duplicate of curved edge "
                    << edges_[cei] << endl;
                ok = false;
                break;
            }
        }
    }

    // Check curved-edge/block-edge correspondence
    //
    // Loop over the edges of each block rather than the edgeList of the
    // topological mesh due to problems with calcEdges for blocks with
    // repeated point labels
    const blockList& blocks = *this;

    for (const blockEdge& curvEdge : edges_)
    {
        bool found = false;

        for (const block& blk : blocks)
        {
            for (const edge& blkEdge : blk.blockShape().edges())
            {
                found = curvEdge.compare(blkEdge) != 0;
                if (found) break;
            }
            if (found) break;
        }

        if (!found)
        {
            Info<< "    Curved edge ";
            curvEdge.write(Info, dict);
            Info<< "    does not correspond to a block edge." << endl;
            ok = false;
        }
    }

    const faceList& faces = bm.faces();

    // Check for duplicate curved face definitions
    forAll(faces_, cfi)
    {
        for (label cfj=cfi+1; cfj<faces_.size(); cfj++)
        {
            if (faces_[cfi].compare(faces_[cfj]) != 0)
            {
                Info<< "    Curved face ";
                faces_[cfj].write(Info, dict);
                Info<< "    is a duplicate of curved face ";
                faces_[cfi].write(Info, dict);
                Info<< endl;
                ok = false;
                break;
            }
        }
    }

    // Check curved-face/block-face correspondence
    for (const blockFace& cface : faces_)
    {
        const face& cf = cface.vertices();

        bool found = false;

        if (cf.size() == 2)
        {
            const label bi = cf[0];
            const label fi = cf[1];

            found =
            (
                bi >= 0 && bi < blocks.size()
             && fi >= 0 && fi < blocks[bi].blockShape().nFaces()
            );
        }

        if (!found)
        {
            for (const face& bf : faces)
            {
                found = cface.compare(bf) != 0;
                if (found) break;
            }
        }

        if (!found)
        {
            Info<< "    Curved face ";
            cface.write(Info, dict);
            Info<< "    does not correspond to a block face." << endl;
            ok = false;
        }
    }


    const pointField& points = bm.points();
    const cellList& cells = bm.cells();
    const polyPatchList& patches = bm.boundaryMesh();

    label nBoundaryFaces = 0;
    for (const cell& c : cells)
    {
        nBoundaryFaces += c.nFaces();
    }
    nBoundaryFaces -= 2*bm.nInternalFaces();

    label nDefinedBoundaryFaces = 0;
    for (const polyPatch& pp : patches)
    {
        nDefinedBoundaryFaces += pp.size();
    }


    if (verbose_)
    {
        Info<< nl << tab << "Basic statistics" << nl
            << tab << tab << "Number of internal faces : "
            << bm.nInternalFaces() << nl
            << tab << tab << "Number of boundary faces : "
            << nBoundaryFaces << nl
            << tab << tab << "Number of defined boundary faces : "
            << nDefinedBoundaryFaces << nl
            << tab << tab << "Number of undefined boundary faces : "
            << nBoundaryFaces - nDefinedBoundaryFaces << nl;

        if ((nBoundaryFaces - nDefinedBoundaryFaces) > 0)
        {
            Info<< tab << tab << tab
                << "(Warning : only leave undefined the front and back planes "
                << "of 2D planar geometries!)" << endl;
        }

        Info<< tab << "Checking patch -> block consistency" << endl;
    }


    for (const polyPatch& pp : patches)
    {
        forAll(pp, patchFacei)
        {
            const face& patchFace = pp[patchFacei];

            bool patchFaceOK = false;

            for (const labelList& cellFaces : cells)
            {
                for (const label cellFacei : cellFaces)
                {
                    const face& cellFace = faces[cellFacei];

                    if (patchFace == cellFace)
                    {
                        patchFaceOK = true;

                        if
                        (
                            (
                                patchFace.areaNormal(points)
                              & cellFace.areaNormal(points)
                            ) < 0
                        )
                        {
                            Info<< tab << tab
                                << "Face " << patchFacei
                                << " of patch " << pp.index()
                                << " (" << pp.name() << ')'
                                << " points inwards"
                                << endl;

                            ok = false;
                        }
                    }
                }
            }

            if (!patchFaceOK)
            {
                Info<< tab << tab
                    << "Face " << patchFacei
                    << " of patch " << pp.index()
                    << " (" << pp.name() << ')'
                    << " does not match any block faces" << endl;

                ok = false;
            }
        }
    }

    if (verbose_)
    {
        Info<< endl;
    }

    // Report patch block/face correspondence
    if (verbose_ > 1)
    {
        const labelList& own = bm.faceOwner();

        Info.stream().setf(ios_base::left);

        Info<< setw(20) << "patch" << "block/face" << nl
            << setw(20) << "-----" << "----------" << nl;

        for (const polyPatch& pp : patches)
        {
            Info<< setw(20) << pp.name();

            label meshFacei = pp.start();

            forAll(pp, bfacei)
            {
                const label celli = own[meshFacei];
                const label cellFacei = cells[celli].find(meshFacei);

                if (bfacei) Info<< token::SPACE;

                Info<< token::BEGIN_LIST
                    << celli << ' ' << cellFacei
                    << token::END_LIST;

                ++meshFacei;
            }
            Info<< nl;
        }

        Info<< setw(20) << "-----" << "----------" << nl
            << nl;
    }

    if (!ok)
    {
        FatalErrorInFunction
            << "Block mesh topology incorrect, stopping mesh generation!"
            << exit(FatalError);
    }
}


// ************************************************************************* //

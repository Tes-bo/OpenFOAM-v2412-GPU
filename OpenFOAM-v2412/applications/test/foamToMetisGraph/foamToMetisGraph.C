/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

Description
    Create a metis graph file representation of an OpenFOAM mesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addNote
    (
        "Create a metis graph file representation for an OpenFOAM mesh"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const labelListList& cellCells = mesh.cellCells();

    // No. of Nodes = nCells
    // No. of Edges connecting Nodes = nInternalFaces

    OFstream os(args.caseName() + ".graph", IOstreamOption::ASCII);

    os  << "%% metis graph file, of an OpenFOAM mesh %%" << nl
        << "%% nCells=" << mesh.nCells()
        << " nFaces=" << mesh.nFaces()
        << " nInternalFaces=" << mesh.nInternalFaces() << nl;

    os  << cellCells.size() << " " << mesh.nInternalFaces() << nl;

    for (const auto& edges : cellCells)
    {
        forAll(edges, i)
        {
            if (i) os << " ";
            os << edges[i] + 1;  // index starts at 1.
        }
        os << nl;
    }

    Info<<"Wrote graph with "
        << mesh.nCells() << " nodes and "
        << mesh.nInternalFaces() << " edges to "
        << os.name() << nl;

    Info<< nl << "End\n" << endl;

    return 0;
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    surfaceCoarsen

Group
    grpSurfaceUtilities

Description
    Surface coarsening using 'bunnylod'

    Reference:
    \verbatim
        Polygon Reduction Demo
        By Stan Melax (c) 1998
        mailto:melax@cs.ualberta.ca
        http://www.cs.ualberta.ca/~melax
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "triFace.H"
#include "triFaceList.H"

// From bunnylod
#include "progmesh.hxx"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int mapVertex(::List<int>& collapse_map, int a, int mx)
{
    if (mx <= 0)
    {
        return 0;
    }
    while (a >= mx)
    {
        a = collapse_map[a];
    }
    return a;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Surface coarsening using 'bunnylod'"
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("factor", "The reduction factor [0,1)");
    argList::addArgument("output", "The output surface file");
    argList::addOption
    (
        "scale",
        "factor",
        "Input geometry scaling factor"
    );

    argList args(argc, argv);

    const auto inFileName = args.get<fileName>(1);
    const auto reduction = args.get<scalar>(2);
    const auto outFileName = args.get<fileName>(3);

    if (reduction <= 0 || reduction > 1)
    {
        FatalErrorInFunction
            << "Reduction factor " << reduction
            << " should be within 0..1" << endl
            << "(it is the reduction in number of vertices)"
            << exit(FatalError);
    }

    const scalar scaleFactor = args.getOrDefault<scalar>("scale", -1);

    Info<< "Input surface   :" << inFileName << nl
        << "Scaling factor  :" << scaleFactor << nl
        << "Reduction factor:" << reduction << nl
        << "Output surface  :" << outFileName << nl
        << endl;

    const triSurface surf(inFileName, scaleFactor);

    Info<< "Surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    ::List<::Vector> vert;     // global list of vertices
    ::List<::tridata> tri;     // global list of triangles


    // Convert triSurface to progmesh format. Note: can use global point
    // numbering since surface read in from file.
    const pointField& pts = surf.points();

    for (const point& pt : pts)
    {
        vert.Add(::Vector(pt.x(), pt.y(), pt.z()));
    }

    for (const labelledTri& f : surf)
    {
        tridata td;
        td.v[0] = f[0];
        td.v[1] = f[1];
        td.v[2] = f[2];
        tri.Add(td);
    }

    ::List<int> collapse_map;   // to which neighbor each vertex collapses
    ::List<int> permutation;

    ::ProgressiveMesh(vert,tri,collapse_map,permutation);

    // rearrange the vertex list
    ::List<::Vector> temp_list;
    for (int i=0; i<vert.num; i++)
    {
        temp_list.Add(vert[i]);
    }
    for (int i=0; i<vert.num; i++)
    {
        vert[permutation[i]] = temp_list[i];
    }

    // update the changes in the entries in the triangle list
    for (int i=0; i<tri.num; i++)
    {
        for (int j=0; j<3; j++)
        {
            tri[i].v[j] = permutation[tri[i].v[j]];
        }
    }

    // Only get triangles with non-collapsed edges.
    int render_num = int(reduction * surf.nPoints());

    Info<< "Reducing to " << render_num << " vertices" << endl;


    // Storage for new surface.
    Foam::List<labelledTri> newTris(surf.size());

    label newI = 0;

    for (int i=0; i<tri.num; i++)
    {
        int p0 = mapVertex(collapse_map, tri[i].v[0], render_num);
        int p1 = mapVertex(collapse_map, tri[i].v[1], render_num);
        int p2 = mapVertex(collapse_map, tri[i].v[2], render_num);

        // note:  serious optimization opportunity here,
        //  by sorting the triangles the following "continue"
        //  could have been made into a "break" statement.
        if (p0 == p1 || p1 == p2 || p2 == p0)
        {
            continue;
        }

        newTris[newI++] = labelledTri(p0, p1, p2, 0);
    }
    newTris.setSize(newI);

    // Convert vert into pointField.
    pointField newPoints(vert.num);

    for (int i=0; i<vert.num; i++)
    {
        const ::Vector & v = vert[i];

        newPoints[i] = point(v.x, v.y, v.z);
    }

    triSurface surf2(newTris, newPoints);

    triSurface outSurf
    (
        surf2.localFaces(),
        surf2.patches(),
        surf2.localPoints()
    );

    Info<< "Coarsened surface:" << endl;
    surf2.writeStats(Info);
    Info<< endl;

    Info<< "Writing to file " << outFileName << endl << endl;

    surf2.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

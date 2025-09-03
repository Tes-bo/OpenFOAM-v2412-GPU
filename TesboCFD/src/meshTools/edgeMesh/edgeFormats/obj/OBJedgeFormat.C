/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "OBJedgeFormat.H"
#include "clock.H"
#include "Fstream.H"
#include "Ostream.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Token list with one of the following:
//     f v1 v2 v3 ...
//     f v1/vt1 v2/vt2 v3/vt3 ...
//     l v1 v2 v3 ...
//     l v1/vt1 v2/vt2 v3/vt3 ...
static label readObjVertices
(
    const SubStrings<string>& tokens,
    DynamicList<label>& verts
)
{
    verts.clear();

    bool first = true;
    for (const auto& tok : tokens)
    {
        if (first)
        {
            // skip initial "f" or "l"
            first = false;
            continue;
        }

        const string vrtSpec(tok);
        const auto slash = vrtSpec.find('/');

        const label vertId =
        (
            slash != string::npos
          ? readLabel(vrtSpec.substr(0, slash))
          : readLabel(vrtSpec)
        );

        verts.append(vertId - 1);
    }

    return verts.size();
}

} // End namespace Foam



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::OBJedgeFormat::OBJedgeFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileFormats::OBJedgeFormat::read(const fileName& filename)
{
    clear();

    IFstream is(filename);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << filename
            << exit(FatalError);
    }

    DynamicList<point> dynPoints;
    DynamicList<label> dynVerts;
    DynamicList<edge>  dynEdges;
    DynamicList<label> dynUsedPoints;

    while (is.good())
    {
        string line = this->getLineNoComment(is);

        // Line continuations
        while (line.removeEnd('\\'))
        {
            line += this->getLineNoComment(is);
        }

        const auto tokens = stringOps::splitSpace(line);

        // Require command and some arguments
        if (tokens.size() < 2)
        {
            continue;
        }

        const word cmd = word::validate(tokens[0]);

        if (cmd == "v")
        {
            // Vertex
            // v x y z

            dynPoints.append
            (
                point
                (
                    readScalar(tokens[1]),
                    readScalar(tokens[2]),
                    readScalar(tokens[3])
                )
            );

            dynUsedPoints.append(-1);
        }
        else if (cmd == "l")
        {
            // Line
            // l v1 v2 v3 ...
            // OR
            // l v1/vt1 v2/vt2 v3/vt3 ...

            readObjVertices(tokens, dynVerts);

            for (label i = 1; i < dynVerts.size(); i++)
            {
                const edge e(dynVerts[i-1], dynVerts[i]);
                dynEdges.append(e);

                dynUsedPoints[e[0]] = e[0];
                dynUsedPoints[e[1]] = e[1];
            }
        }
        else if (cmd == "f")
        {
            // Face - support for faces with 2 vertices

            // f v1 v2 v3 ...
            // OR
            // f v1/vt1 v2/vt2 v3/vt3 ...

            if (readObjVertices(tokens, dynVerts) == 2)
            {
                for (label i = 1; i < dynVerts.size(); i++)
                {
                    const edge e(dynVerts[i-1], dynVerts[i]);
                    dynEdges.append(e);

                    dynUsedPoints[e[0]] = e[0];
                    dynUsedPoints[e[1]] = e[1];
                }
            }
        }
    }

    // Cull unused points
    label nUsed = 0;
    forAll(dynPoints, pointi)
    {
        if (dynUsedPoints[pointi] >= 0)
        {
            if (nUsed != pointi)
            {
                dynPoints[nUsed] = std::move(dynPoints[pointi]);
                dynUsedPoints[pointi] = nUsed;   // The new list location
            }
            ++nUsed;
        }
    }

    dynPoints.setSize(nUsed);

    // Transfer to normal lists
    storedPoints().transfer(dynPoints);

    // Renumber edge vertices
    if (nUsed != dynUsedPoints.size())
    {
        for (edge& e : dynEdges)
        {
            e[0] = dynUsedPoints[e[0]];
            e[1] = dynUsedPoints[e[1]];
        }
    }
    storedEdges().transfer(dynEdges);

    return true;
}


void Foam::fileFormats::OBJedgeFormat::write
(
    const fileName& filename,
    const edgeMesh& mesh,
    IOstreamOption streamOpt,
    const dictionary&
)
{
    // ASCII only, allow output compression
    streamOpt.format(IOstreamOption::ASCII);

    const pointField& pointLst = mesh.points();
    const edgeList& edgeLst = mesh.edges();

    OFstream os(filename, streamOpt);
    if (!os.good())
    {
        FatalErrorInFunction
            << "Cannot write file " << filename << nl
            << exit(FatalError);
    }


    os  << "# Wavefront OBJ file written " << clock::dateTime().c_str() << nl
        << "o " << os.name().stem() << nl
        << nl
        << "# points : " << pointLst.size() << nl
        << "# lines  : " << edgeLst.size() << nl;

    os  << nl
        << "# <points count=\"" << pointLst.size() << "\">" << nl;

    // Write vertex coords
    for (const point& p : pointLst)
    {
        os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
    }

    os  << "# </points>" << nl
        << nl
        << "# <edges count=\"" << edgeLst.size() << "\">" << endl;

    // Write line connectivity
    for (const edge& e : edgeLst)
    {
        os << "l " << (e[0] + 1) << " " << (e[1] + 1) << nl;
    }
    os << "# </edges>" << endl;
}


// ************************************************************************* //

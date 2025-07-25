/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    objToVTK

Group
    grpMeshManipulationUtilities

Description
    Read obj line (not surface) file and convert into legacy VTK file.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OFstream.H"
#include "StringStream.H"
#include "point.H"
#include "DynamicList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string getLine(std::ifstream& is)
{
    string line;
    do
    {
        std::getline(is, line);
    }
    while (line.starts_with('#'));

    return line;
}


// Read space-separated vertices (with optional '/' arguments)
labelList parseVertices(const string& line)
{
    DynamicList<label> verts;

    // Assume 'l' is followed by space.
    string::size_type endNum = 1;

    do
    {
        string::size_type startNum = line.find_first_not_of(' ', endNum);

        if (startNum == string::npos)
        {
            break;
        }

        endNum = line.find(' ', startNum);

        string vertexSpec;
        if (endNum != string::npos)
        {
            vertexSpec = line.substr(startNum, endNum-startNum);
        }
        else
        {
            vertexSpec = line.substr(startNum, line.size() - startNum);
        }

        string::size_type slashPos = vertexSpec.find('/');

        label vertI = 0;
        if (slashPos != string::npos)
        {
            IStringStream intStream(vertexSpec.substr(0, slashPos));

            intStream >> vertI;
        }
        else
        {
            IStringStream intStream(vertexSpec);

            intStream >> vertI;
        }
        verts.append(vertI - 1);
    }
    while (true);

    return verts.shrink();
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Read obj line (not surface) file and convert into legacy VTK file"
    );

    argList::noParallel();
    argList::addArgument("obj-file", "The input obj line file");
    argList::addArgument("vtk-file", "The output vtk file");
    argList args(argc, argv);

    const auto objName = args.get<fileName>(1);
    const auto outName = args.get<fileName>(2);

    std::ifstream OBJfile(objName);

    Info<< "Processing file " << objName << endl;

    if (!OBJfile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << objName << exit(FatalError);
    }

    // Points and lines
    DynamicList<point> points;
    DynamicList<vector> pointNormals;
    DynamicList<labelList> polyLines;
    DynamicList<labelList> polygons;

    bool hasWarned = false;

    label lineNo = 0;
    while (OBJfile.good())
    {
        const string line = getLine(OBJfile);
        lineNo++;

        if (line.empty()) continue;

        // Read first word
        IStringStream lineStream(line);
        word cmd(lineStream);

        if (cmd == "v")
        {
            scalar x, y, z;

            lineStream >> x >> y >> z;

            points.append(point(x, y, z));
        }
        else if (cmd == "vn")
        {
            scalar x, y, z;

            lineStream >> x >> y >> z;

            pointNormals.append(vector(x, y, z));
        }
        else if (cmd == "l")
        {
            polyLines.append(parseVertices(line));
        }
        else if (cmd == "f")
        {
            polygons.append(parseVertices(line));
        }
        else if (cmd != "")
        {
            if (!hasWarned)
            {
                hasWarned = true;

                WarningInFunction
                    << "Unrecognized OBJ command " << cmd << nl
                    << "In line " << lineStream.str()
                    << " at linenumber " << lineNo << nl
                    << "Only recognized commands are 'v' and 'l'.\n"
                    << "If this is a surface command use surfaceConvert instead"
                    << " to convert to a file format that can be read by VTK"
                    << endl;
            }
        }
    }


    //
    // Write as vtk 'polydata' file
    //


    OFstream outFile(outName);

    outFile
        << "# vtk DataFile Version 2.0\n"
        << objName << nl
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << points.size() << " double\n";

    for (const point& pt : points)
    {
        outFile
            << float(pt.x()) << ' '
            << float(pt.y()) << ' '
            << float(pt.z()) << nl;
    }

    outFile
        << "VERTICES " << points.size() << ' ' << (2 * points.size()) << nl;

    forAll(points, i)
    {
        outFile << 1 << ' ' << i << nl;
    }

    label nItems = 0;
    forAll(polyLines, polyI)
    {
        nItems += polyLines[polyI].size() + 1;
    }

    outFile
        << "LINES " << polyLines.size() << ' ' << nItems << nl;

    forAll(polyLines, polyI)
    {
        const labelList& line = polyLines[polyI];

        outFile << line.size();

        forAll(line, i)
        {
            outFile << ' ' << line[i];
        }
        outFile << nl;
    }


    nItems = 0;
    forAll(polygons, polyI)
    {
        nItems += polygons[polyI].size() + 1;
    }

    outFile
        << "POLYGONS " << polygons.size() << ' ' << nItems << nl;

    forAll(polygons, polyI)
    {
        const labelList& line = polygons[polyI];

        outFile << line.size();

        forAll(line, i)
        {
            outFile << ' ' << line[i];
        }
        outFile << nl;
    }


    outFile
        << "POINT_DATA " << points.size() << nl
        << "SCALARS pointID double 1\n"
        << "LOOKUP_TABLE default\n";

    forAll(points, i)
    {
        outFile << i;

        if ((i % 10) == 1)
        {
            outFile << nl;
        }
        else
        {
            outFile << ' ';
        }
    }

    if (!pointNormals.empty())
    {
        outFile << nl << "NORMALS pointNormals double\n";

        for(const vector& n : pointNormals)
        {
            outFile
                << float(n.x()) << ' '
                << float(n.y()) << ' '
                << float(n.z()) << nl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

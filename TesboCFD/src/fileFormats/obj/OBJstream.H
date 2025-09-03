/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Class
    Foam::OBJstream

Description
    An OFstream that keeps track of vertices and provides convenience
    output methods for OBJ files.

SourceFiles
    OBJstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_OBJstream_H
#define Foam_OBJstream_H

#include "OFstream.H"
#include "point.H"
#include "edge.H"
#include "face.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class treeBoundBox;

/*---------------------------------------------------------------------------*\
                          Class OBJstream Declaration
\*---------------------------------------------------------------------------*/

class OBJstream
:
    public OFstream
{
    // Private Data

        bool startOfLine_;

        label nVertices_;


    // Private Member Functions

        //- State engine to count number of vertices.
        //- Triggered on newline and 'v' (vertex).
        inline void vertex_state(const char c);

public:

    //- Declare type-name (no debug switch)
    ClassNameNoDebug("OBJstream");


    // Constructors

        //- Construct from pathname (ASCII, uncompressed)
        explicit OBJstream
        (
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct from pathname, format (uncompressed)
        OBJstream
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            OBJstream(pathname, IOstreamOption(fmt, cmp))
        {}


    //- Destructor
    ~OBJstream() = default;


    // Member Functions

    // Access

        //- Return the number of vertices written
        label nVertices() const noexcept { return nVertices_; }


    // Write Functions

        //- Inherit write methods from OFstream
        using OFstream::write;
        using OFstream::writeQuoted;

        //- Write character
        virtual Ostream& write(const char c) override;

        //- Write character/string content, with/without surrounding quotes
        virtual Ostream& writeQuoted
        (
            const char* str,
            std::streamsize len,
            const bool quoted=true
        ) override;

        //- Write character string
        virtual Ostream& write(const char* str) override;

        //- Write word
        virtual Ostream& write(const word& str) override;

        //- Write string
        virtual Ostream& write(const std::string& str) override;


    // Direct write functionality

        //- Write comment (with '# ' prefix)
        Ostream& writeComment(const std::string& str);

        //- Write point (vertex)
        Ostream& write(const point& p);

        //- Write point and vector normal ('vn')
        Ostream& write(const point& p, const vector& n);

        //- Write multiple points
        Ostream& write(const UList<point>& points);

        //- Write edge as points with line
        Ostream& write(const edge& e, const UList<point>& points);

        //- Write line
        Ostream& write(const linePointRef& ln);

        //- Write line with points and vector normals ('vn')
        Ostream& write
        (
            const linePointRef& ln,
            const vector& n0,
            const vector& n1
        );

        //- Write line joining two points
        Ostream& writeLine(const point& p0, const point& p1);

        //- Write triangle as points with lines/filled-polygon
        Ostream& write(const triPointRef& f, const bool lines = true);

        //- Write face loop points with lines/filled-polygon
        Ostream& writeFace
        (
            const UList<point>& points,
            const bool lines = true
        );

        //- Write face as points with lines/filled-polygon
        Ostream& write
        (
            const face& f,
            const UList<point>& points,
            const bool lines = true
        );

        //- Write patch faces as points with lines/filled-polygon
        Ostream& write
        (
            const UList<face>& faces,
            const pointField& points,
            const bool lines = true
        );

        //- Write edges as points with lines.
        //  Optionally eliminate unused points.
        Ostream& write
        (
            const UList<edge>& edges,
            const UList<point>& points,
            const bool compact = false
        );

        //- Write tree-bounding box with lines/filled-polygons
        Ostream& write(const treeBoundBox& bb, const bool lines = true);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

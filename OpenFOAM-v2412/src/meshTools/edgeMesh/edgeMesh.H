/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::edgeMesh

Description
    Points connected by edges.

    Can be read from fileName based on extension. Uses ::New factory method
    to select the reader and transfer the result.

SourceFiles
    edgeMeshI.H
    edgeMesh.C
    edgeMeshIO.C
    edgeMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_edgeMesh_H
#define Foam_edgeMesh_H

#include "pointField.H"
#include "edgeList.H"
#include "dictionary.H"
#include "edgeMeshFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class edgeMesh;

Istream& operator>>(Istream& is, edgeMesh& em);
Ostream& operator<<(Ostream& os, const edgeMesh& em);

/*---------------------------------------------------------------------------*\
                           Class edgeMesh Declaration
\*---------------------------------------------------------------------------*/

class edgeMesh
:
    public fileFormats::edgeMeshFormatsCore
{
    // Private Data

        //- Vertices of the edges
        pointField points_;

        //- The edges defining the boundary
        edgeList edges_;

        //- From point to edges
        mutable std::unique_ptr<labelListList> pointEdgesPtr_;


    // Private Member Functions

        //- Calculate point-edge addressing (inverse of edges)
        void calcPointEdges() const;


protected:

    // Protected Member Functions

        //- Non-const access to global points
        inline pointField& storedPoints() noexcept;

        //- Non-const access to the edges
        inline edgeList& storedEdges() noexcept;


public:

    //- Runtime type information
    TypeName("edgeMesh");


    // Static Member Functions

        //- Summary of supported read file types.
        static wordHashSet readTypes();

        //- Summary of supported write file types.
        static wordHashSet writeTypes();

        //- Can we read this file format?
        static bool canReadType(const word& fileType, bool verbose=false);

        //- Can we write this file format type?
        static bool canWriteType(const word& fileType, bool verbose=false);

        //- Can we read this file format?
        static bool canRead(const fileName& name, bool verbose=false);


    // Constructors

        //- Default construct
        inline edgeMesh();

        //- Copy construct
        inline edgeMesh(const edgeMesh& em);

        //- Move construct
        inline edgeMesh(edgeMesh&& em);

        //- Copy construct from components
        inline edgeMesh(const pointField& points, const edgeList& edges);

        //- Move construct from components
        inline edgeMesh(pointField&& pointLst, edgeList&& edgeLst);

        //- Construct from file name (uses extension to determine type)
        explicit edgeMesh(const fileName& name);

        //- Construct from file name with specified type
        edgeMesh(const fileName& name, const word& fileType);


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            edgeMesh,
            fileExtension,
            (
                const fileName& name
            ),
            (name)
        );


    // Selectors

        //- Read construct from filename with given format
        static autoPtr<edgeMesh> New
        (
            const fileName& name,
            const word& fileType
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<edgeMesh> New(const fileName& name);


    //- Destructor
    virtual ~edgeMesh() = default;


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            edgeMesh,
            write,
            fileExtension,
            (
                const fileName& name,
                const edgeMesh& mesh,
                IOstreamOption streamOpt,
                const dictionary& options
            ),
            (name, mesh, streamOpt, options)
        );

        //- Write to file (format implicit in the extension)
        static void write
        (
            const fileName& name,
            const edgeMesh& mesh,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );

        //- Write to file, with given format
        static void write
        (
            const fileName& name,
            const word& fileType,
            const edgeMesh& mesh,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        );


    // Member Functions

        //- Transfer the contents of the argument and annul the argument
        void transfer(edgeMesh& mesh);


    // Read

        //- Read from file. Chooses reader based on explicit extension
        bool read(const fileName& name, const word& fileType);

        //- Read from file. Chooses reader based on detected extension
        virtual bool read(const fileName& name);


    // Access

        //- Return points
        inline const pointField& points() const noexcept;

        //- Return edges
        inline const edgeList& edges() const noexcept;

        //- Return edges
        inline const labelListList& pointEdges() const;

        //- Find connected regions. Set region number per edge.
        //  Returns number of regions.
        label regions(labelList& edgeRegion) const;


    // Edit

        //- Clear all storage
        virtual void clear();

        //- Scale points. A non-positive factor is ignored
        virtual void scalePoints(const scalar scaleFactor);

        //- Geometric merge points (points within mergeDist) prior to
        //  automatically calling mergeEdges().
        virtual void mergePoints(const scalar mergeDist);

        //- Merge duplicate edges and eliminate unused points.
        virtual void mergeEdges();


    // Write

        virtual void writeStats(Ostream&) const;

        //- Write to file, choosing writer based on the file extension.
        virtual void write
        (
            const fileName& name,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, *this, streamOpt, options);
        }

        //- Write to file with given format type.
        //  If the format type is "", uses the file extension.
        virtual void write
        (
            const fileName& name,
            const word& fileType,
            IOstreamOption streamOpt = IOstreamOption(),
            const dictionary& options = dictionary::null
        ) const
        {
            write(name, fileType, *this, streamOpt, options);
        }


    // Member Operators

        //- Copy assignment
        inline void operator=(const edgeMesh& rhs);

        //- Move assignment
        inline void operator=(edgeMesh&& rhs);


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const edgeMesh& em);
        friend Istream& operator>>(Istream& is, edgeMesh& em);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

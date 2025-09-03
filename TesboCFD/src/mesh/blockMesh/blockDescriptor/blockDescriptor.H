/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::blockDescriptor

Description
    Takes the description of the block and the list of curved edges and
    creates a list of points on edges together with the weighting factors

    For a given block, the correspondence between the ordering of vertex labels
    and face labels is shown below.  For vertex numbering in the sequence 0 to 7
    (block, centre): faces 0 (f0) and 1 are left and right, respectively; faces
    2 and 3 are front and back; and faces 4 and 5 are bottom and top:

    \verbatim
                    7 ---- 6
              f5    |\     :\    f3
               |    | 4 ---- 5    \
               |    3.|....2 |     \
               |     \|     \|     f2
              f4      0 ---- 1
      Y  Z
       \ |          f0 ------ f1
        \|
         o--- X
     \endverbatim

SourceFiles
    blockDescriptor.C
    blockDescriptorEdges.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_blockDescriptor_H
#define Foam_blockDescriptor_H

#include "ijkMesh.H"
#include "cellShape.H"
#include "pointField.H"
#include "scalarList.H"
#include "blockEdgeList.H"
#include "blockFaceList.H"
#include "gradingDescriptors.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class blockDescriptor;

Ostream& operator<<(Ostream&, const blockDescriptor&);

Ostream& operator<<(Ostream&, const InfoProxy<blockDescriptor>&);


/*---------------------------------------------------------------------------*\
                       Class blockDescriptor Declaration
\*---------------------------------------------------------------------------*/

class blockDescriptor
:
    public ijkMesh
{
    // Private Data

        //- Reference to point field defining the block mesh
        const pointField& vertices_;

        //- Reference to a list of block edges
        const blockEdgeList& blockEdges_;

        //- Reference to the list of curved faces
        const blockFaceList& blockFaces_;

        //- Block shape. Likely only hex is supportable
        cellShape blockShape_;

        //- Expansion ratios in all directions
        List<gradingDescriptors> expand_;

        //- The block index in the originating list (-1 for unknown)
        label index_;

        //- Name of the zone (empty word if none)
        word zoneName_;

        //- Curved-face labels for each block-face (-1 for flat faces)
        FixedList<label, 6> curvedFaces_;

        //- Number of curved faces in this block
        label nCurvedFaces_;


    // Private Member Functions

        //- Assign edge grading.
        //  \return false for unsupported specification
        bool assignGradings(const UList<gradingDescriptors>& ratios);

        //- Check block has outward-pointing faces
        void check(const Istream& is);

        //- Calculate the points and weights for the specified hex edge.
        //  Return the number of curved edges (0-1)
        int calcEdgePointsWeights
        (
            pointField& edgePoints,
            scalarList& edgeWeights,
            const Foam::edge& cellModelEdge,
            const label nDiv,
            const gradingDescriptors& expand
        ) const;

        void findCurvedFaces(const label blockIndex = -1);

public:

    // Generated Methods

        //- Copy construct
        blockDescriptor(const blockDescriptor&) = default;

        //- No copy assignment
        void operator=(const blockDescriptor&) = delete;


    // Constructors

        //- Construct from components. Optional zone name.
        blockDescriptor
        (
            const cellShape& bshape,
            const pointField& vertices,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            const labelVector& density,
            const UList<gradingDescriptors>& expand
                = UList<gradingDescriptors>::null(),
            const word& zoneName = word::null
        );

        //- Construct from Istream
        blockDescriptor
        (
            const dictionary& dict,
            const label blockIndex,
            const pointField& vertices,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            Istream& is
        );


    // Member Functions

        //- Reference to point field defining the block mesh
        inline const pointField& vertices() const noexcept;

        //- Return reference to the list of curved faces
        inline const blockFaceList& blockFaces() const noexcept;

        //- Return the block shape
        inline const cellShape& blockShape() const noexcept;

        //- The mesh density (number of cells) in the i,j,k directions
        inline const labelVector& density() const noexcept;

        //- Expansion ratios in all directions
        inline const List<gradingDescriptors>& grading() const noexcept;

        //- Return the (optional) zone name
        inline const word& zoneName() const noexcept;

        //- Curved-face labels for each block-face (-1 for flat faces)
        inline const FixedList<label, 6>& curvedFaces() const noexcept;

        //- Number of curved faces in this block
        inline label nCurvedFaces() const noexcept;

        //- Return block point for local label i
        inline const point& blockPoint(const label i) const;

        //- Face vertex label offset for a particular i,j,k position
        //- on hex face (0-5)
        inline label facePointLabel
        (
            const direction facei,
            const label i,
            const label j
        ) const;

        //- True if point i,j,k addresses a block vertex
        inline bool vertex(const label i, const label j, const label k) const;

        //- True if point i,j,k addresses a block edge
        inline bool edge(const label i, const label j, const label k) const;

        //- Calculate the points and weights for all edges.
        //  \return the number of curved edges (0-12)
        int edgesPointsWeights
        (
            pointField (&edgesPoints)[12],
            scalarList (&edgesWeights)[12]
        ) const;

        //- Calculate points and weights for specified edge,
        //- using the specified number of divisions and grading
        //  \return True if the edge is curved
        bool edgePointsWeights
        (
            const label edgei,
            pointField& edgePoints,
            scalarList& edgeWeights,
            const label nDiv,
            const gradingDescriptors& gd = gradingDescriptors()
        ) const;

        //- Calculate points and weights for specified edge.
        //  \return True if the edge is curved
        bool edgePointsWeights
        (
            const label edgei,
            pointField& edgePoints,
            scalarList& edgeWeights
        ) const;

        //- Return true if point i,j,k addresses a block flat face or edge
        inline bool flatFaceOrEdge
        (
            const label i,
            const label j,
            const label k
        ) const;

        //- Return the list of face-points for all of the faces of the block
        FixedList<pointField, 6> facePoints(const pointField& points) const;

        //- Correct the location of the given face-points
        //  to lie on the faces of the block
        void correctFacePoints(FixedList<pointField, 6>&) const;

        //- Write block index with dictionary lookup
        static void write(Ostream&, const label blocki, const dictionary&);

        //- Return info proxy,
        //- used to print information to a stream
        InfoProxy<blockDescriptor> info() const noexcept { return *this; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blockDescriptorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::block

Description
    Creates a single block of cells from point coordinates, numbers of
    cells in each direction and an expansion ratio.

Note
    The cells for filling the block are demand-driven.

SourceFiles
    block.C
    blockCreate.C

\*---------------------------------------------------------------------------*/

#ifndef block_H
#define block_H

#include "blockDescriptor.H"
#include "cellShapeList.H"
#include "hexCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class block Declaration
\*---------------------------------------------------------------------------*/

class block
:
    public blockDescriptor
{
    // Private Data

        //- List of points
        pointField points_;

        //- The cells (hex)
        List<hexCell> blockCells_;

        //- Boundary patches
        FixedList<List<FixedList<label, 4>>, 6> blockPatches_;


    // Private Member Functions

        //- Create vertices for cells filling the block
        void createPoints();

        //- Create cells
        void createCells();

        //- Create boundary patch faces for the block
        void createBoundary();

        //- Add boundary faces for the shape face to the output list at
        //- the iterator location
        template<class OutputIterator>
        OutputIterator addBoundaryFaces
        (
            const direction shapeFacei,
            OutputIterator iter
        ) const;


        //- No copy construct
        block(const block&) = delete;

        //- No copy assignment
        void operator=(const block&) = delete;


public:

    //- Runtime type information
    TypeName("block");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            block,
            Istream,
            (
                const dictionary& dict,
                const label index,
                const pointField& vertices,
                const blockEdgeList& edges,
                const blockFaceList& faces,
                Istream& is
            ),
            (dict, index, vertices, edges, faces, is)
        );


    // Constructors

        //- Construct from components. Optional zone name.
        block
        (
            const cellShape& bshape,
            const pointField& vertices,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            const labelVector& density,
            const UList<gradingDescriptors>& expand
                = UList<gradingDescriptors>::null(),
            const word& zoneName = ""
        );

        //- Construct from components with Istream
        block
        (
            const dictionary& dict,
            const label index,
            const pointField& vertices,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            Istream& is
        );

        //- Construct from a block definition
        block(const blockDescriptor& blockDesc);

        //- Clone
        autoPtr<block> clone() const
        {
            NotImplemented;
            return nullptr;
        }

        //- New function which constructs and returns pointer to a block
        static autoPtr<block> New
        (
            const dictionary& dict,
            const label index,
            const pointField& points,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            Istream&
        );

        //- Class used for the read-construction of
        //  PtrLists of blocks
        class iNew
        {
            const dictionary& dict_;
            const pointField& points_;
            const blockEdgeList& edges_;
            const blockFaceList& faces_;
            mutable label index_;

        public:

            iNew
            (
                const dictionary& dict,
                const pointField& points,
                const blockEdgeList& edges,
                const blockFaceList& faces
            )
            :
                dict_(dict),
                points_(points),
                edges_(edges),
                faces_(faces),
                index_(0)
            {}

            autoPtr<block> operator()(Istream& is) const
            {
                return block::New(dict_, index_++, points_, edges_, faces_, is);
            }
        };


    //- Destructor
    virtual ~block() = default;


    // Member Functions

    // Access

        //- The points for filling the block
        inline const pointField& points() const noexcept;

        //- The hex cells for filling the block
        inline const List<hexCell>& cells() const;

        //- The boundary patch faces for the block
        inline const FixedList<List<FixedList<label, 4>>, 6>&
        boundaryPatches() const noexcept;


    // Mesh Components

        //- The (hex) cell shapes for filling the block.
        cellShapeList shapes() const;
};


//- Ostream Operator
Ostream& operator<<(Ostream& os, const block& blk);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blockI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

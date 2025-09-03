/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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
    Foam::structuredRenumber

Description
    Renumbering according to mesh layers.
    depthFirst = true:
        first column gets ids 0..nLayer-1,
        second nLayers..2*nLayers-1 etc.
    depthFirst = false:
        first layer gets ids 0,1,2 etc.

SourceFiles
    structuredRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_structuredRenumber_H
#define Foam_structuredRenumber_H

#include "renumberMethod.H"

namespace Foam
{

// Forward Declarations
template<class Type> class topoDistanceData;

/*---------------------------------------------------------------------------*\
                     Class structuredRenumber Declaration
\*---------------------------------------------------------------------------*/

class structuredRenumber
:
    public renumberMethod
{
public:

    // Public Classes

        //- Function class for sorting according to column and layer
        class layerLess
        {
            const bool depthFirst_;
            const labelList& order_;
            const List<topoDistanceData<label>>& distance_;

        public:

            layerLess
            (
                const bool depthFirst,
                const labelList& order,
                const List<topoDistanceData<label>>& distance
            )
            :
                depthFirst_(depthFirst),
                order_(order),
                distance_(distance)
            {}

            bool operator()(const label a, const label b) const;
        };


private:

    // Private Data

        const dictionary& coeffsDict_;

        const wordRes patches_;

        const label nLayers_;

        const bool depthFirst_;

        const bool reverse_;

        const autoPtr<renumberMethod> method_;


    // Private Member Functions

        //- No copy construct
        structuredRenumber(const structuredRenumber&) = delete;

        //- No copy assignment
        void operator=(const structuredRenumber&) = delete;


public:

    //- Runtime type information
    TypeName("structured");


    // Constructors

        //- Construct given the renumber dictionary
        explicit structuredRenumber(const dictionary& dict);


    //- Destructor
    virtual ~structuredRenumber() = default;


    // Member Functions

        //- Renumbering method requires a polyMesh!
        virtual bool needs_mesh() const { return true; }


    // With mesh topology

        //- Return the cell visit order (from ordered back to original cell id)
        //  using the mesh.
        virtual labelList renumber
        (
            const polyMesh& mesh
        ) const;


    // With explicit topology - Not implemented!

        //- Return the cell visit order (from ordered back to original cell id).
        //- Not implemented!
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells
        ) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return the cell visit order (from ordered back to original cell id).
        //- Not implemented!
        virtual labelList renumber
        (
            const labelListList& cellCells
        ) const
        {
            NotImplemented;
            return labelList();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

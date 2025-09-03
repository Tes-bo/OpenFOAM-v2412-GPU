/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2022-2024 OpenCFD Ltd.
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
    Foam::springRenumber

Description
    Use spring analogy - attract neighbouring cells according to the distance
    of their cell indices.

    // Maximum jump of cell indices. Is fraction of number of cells
    maxCo 0.1;

    // Limit the amount of movement; the fraction maxCo gets decreased
    // with every iteration.
    freezeFraction 0.9;

    // Maximum number of iterations
    maxIter 1000;

SourceFiles
    springRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_springRenumber_H
#define Foam_springRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class springRenumber Declaration
\*---------------------------------------------------------------------------*/

class springRenumber
:
    public renumberMethod
{
    // Private Data

        const dictionary& coeffsDict_;

        const label maxIter_;

        const scalar maxCo_;

        const scalar freezeFraction_;

        //- Verbose reporting of progress
        bool verbose_;


    // Private Member Functions

        //- The spring renumber implementation.
        //  The connections are CompactListList<label> or a labelListList.
        template<class ConnectionListListType>
        labelList renumberImpl
        (
            const ConnectionListListType& cellCellAddressing
        ) const;

        //- No copy construct
        springRenumber(const springRenumber&) = delete;

        //- No copy assignment
        void operator=(const springRenumber&) = delete;


public:

    //- Runtime type information
    TypeName("spring");


    // Constructors

        //- Construct given the renumber dictionary
        explicit springRenumber(const dictionary& dict);


    //- Destructor
    virtual ~springRenumber() = default;


    // Member Functions

    // With mesh topology

        //- Return the cell visit order (from ordered back to original cell id)
        //- using the mesh to determine the connectivity.
        virtual labelList renumber
        (
            //! Mesh connectivity (see globalMeshData::calcCellCells)
            const polyMesh& mesh
        ) const;


    // With explicit topology

        //- Return the cell visit order (from ordered back to original cell id)
        virtual labelList renumber
        (
            //! Mesh connectivity
            const CompactListList<label>& cellCells
        ) const;

        //- Return the cell visit order (from ordered back to original cell id)
        virtual labelList renumber
        (
            //! Mesh connectivity
            const labelListList& cellCells
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

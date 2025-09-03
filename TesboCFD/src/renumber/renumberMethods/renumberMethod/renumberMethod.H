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
    Foam::renumberMethod

Description
    Abstract base class for renumbering

SourceFiles
    renumberMethod.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_renumberMethod_H
#define Foam_renumberMethod_H

#include "polyMesh.H"
#include "pointField.H"
#include "CompactListList.H"
#include "wordList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class renumberMethod Declaration
\*---------------------------------------------------------------------------*/

class renumberMethod
{
public:

    //- Runtime type information
    TypeNameNoDebug("renumberMethod");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            renumberMethod,
            dictionary,
            (
                const dictionary& renumberDict
            ),
            (renumberDict)
        );


    // Constructors

        //- Default construct
        renumberMethod()
        {}

        //- Construct with renumber dictionary (which is currently unused)
        explicit renumberMethod(const dictionary&)
        {}


    // Selectors

        //- Construct/select a renumbering method
        static autoPtr<renumberMethod> New(const dictionary& dict);

        //- Return a list of the known methods
        static wordList supportedMethods();


    //- Destructor
    virtual ~renumberMethod() = default;


    // Member Functions

        //- Renumbering method without mesh or cell-cell topology
        //- (very special case)
        virtual bool no_topology() const { return false; }

        //- Renumbering method requires a polyMesh for its topology
        virtual bool needs_mesh() const { return false; }


    // No topology

        //- Return the cell visit order (from ordered back to original cell id)
        //- based solely on the number of cells.
        //  Only applicable for no_topology() methods.
        virtual labelList renumber(const label nCells) const;

        //- Return the cell visit order (from ordered back to original cell id)
        //- based solely on the cell centres (or number of cell centres).
        //  Only applicable for no_topology() methods.
        virtual labelList renumber(const pointField&) const;


    // With mesh topology

        //- Return the cell visit order (from ordered back to original cell id)
        //- using the mesh to determine the connectivity.
        virtual labelList renumber(const polyMesh& mesh) const;


    // With explicit topology

        //- Return the cell visit order (from ordered back to original cell id),
        virtual labelList renumber
        (
            //! Mesh connectivity
            const CompactListList<label>& cellCells
        ) const = 0;

        //- Return the cell visit order (from ordered back to original cell id),
        virtual labelList renumber
        (
            //! Mesh connectivity
            const labelListList& cellCells
        ) const = 0;


    // Housekeeping

        //- Deprecated - the pointField is unused
        //  \deprecated(2024-03) the pointField is unused
        FOAM_DEPRECATED_FOR(2024-03, "renumber(const polyMesh&)")
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField&
        ) const
        {
            return renumber(mesh);
        }

        //- Deprecated - the pointField is unused
        //  \deprecated(2024-03) the pointField is unused
        FOAM_DEPRECATED_FOR(2024-03, "renumber(const CompactListList<label>&)")
        virtual labelList renumber
        (
            const CompactListList<label>& cellCells,
            const pointField&
        ) const
        {
            return renumber(cellCells);
        }

        //- Deprecated - the pointField is unused
        //  \deprecated(2024-03) the pointField is unused
        FOAM_DEPRECATED_FOR(2024-03, "renumber(const labelListList&)")
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField&
        ) const
        {
            return renumber(cellCells);
        }

        //- Deprecated - renumbering with agglomeration map.
        //- Calculate globalCellCells directly
        //  \deprecated(2024-03) calculate globalCellCells directly
        FOAM_DEPRECATED_FOR(2024-03, "calcCellCells and renumber separately")
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const labelUList& fineToCoarse,
            const pointField& coarsePoints
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

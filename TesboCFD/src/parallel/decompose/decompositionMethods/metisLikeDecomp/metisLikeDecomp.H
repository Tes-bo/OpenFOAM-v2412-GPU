/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::metisLikeDecomp

Description
    Domain decomposition using METIS-like data structures.

    When run in parallel will collect the entire graph on to the master,
    decompose and send back.

SourceFiles
    metisLikeDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_metisLikeDecomp_H
#define Foam_metisLikeDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class metisLikeDecomp Declaration
\*---------------------------------------------------------------------------*/

class metisLikeDecomp
:
    public decompositionMethod
{
protected:

    // Protected Data

        //- Coefficients for all derived methods
        const dictionary& coeffsDict_;


    // Protected Member Functions

        //- Serial and/or collect/distribute for parallel operation
        virtual label decomposeGeneral
        (
            const labelList& adjncy,
            const labelList& xadj,
            const List<scalar>& cellWeights,
            labelList& decomp
        ) const;

        //- Decomposition with metis-like parameters
        virtual label decomposeSerial
        (
            const labelList& adjncy,
            const labelList& xadj,
            const List<scalar>& cellWeights,
            labelList& decomp
        ) const = 0;


        //- No copy construct
        metisLikeDecomp(const metisLikeDecomp&) = delete;

        //- No copy assignment
        void operator=(const metisLikeDecomp&) = delete;


    // Constructors

        //- Construct with number of domains (no coefficients or constraints)
        explicit metisLikeDecomp(const label numDomains);

        //- Construct for derived type name and decomposition dictionary.
        metisLikeDecomp
        (
            const word& derivedType,
            const dictionary& decompDict,
            int select = selectionType::NULL_DICT
        )
        :
            metisLikeDecomp(derivedType, decompDict, "", select)
        {}

        //- Construct for derived type name, decomposition dictionary
        //- and (optional) region name
        //  The default search for the coefficients will return dictionary::null
        //  on failure. This avoids a name clash of a metis "method" with the
        //  top level.
        metisLikeDecomp
        (
            const word& derivedType,
            const dictionary& decompDict,
            const word& regionName,
            int select = selectionType::NULL_DICT
        );

public:

    //- Destructor
    virtual ~metisLikeDecomp() = default;


    // Member Functions

        //- Inherit decompose from decompositionMethod
        using decompositionMethod::decompose;

        //- Return for every coordinate the wanted processor number.
        //  Uses the mesh connectivity (if needed).
        //  Weights get normalised so the minimum value is 1 before truncation
        //  to an integer so the weights should be multiples of the minimum
        //  value. The overall sum of weights might otherwise overflow.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points = pointField::null(),
            const scalarField& pointWeights = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Gets passed agglomeration map (from fine to coarse cells) and coarse
        //  cell location. Can be overridden by decomposers that provide this
        //  functionality natively.
        //  See note on weights above.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const labelList& agglom,
            const pointField& agglomPoints,
            const scalarField& agglomWeights = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        //  See note on weights above.
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,
            const pointField& cellCentres = pointField::null(),
            const scalarField& cellWeights = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        //  See note on weights above.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cellCentres = pointField::null(),
            const scalarField& cellWeights = scalarField::null()
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

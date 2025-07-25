/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::randomDecomp

Description
    Decomposition according to pseudo-random number generator,
    which is useful for development and stress testing implementations
    but obviously sub-optimal for normal uses.
    Can optionally specify the number of consecutive cells to assign
    to the same processor (agglomeration).

    Coefficients:
    \table
        Property  | Description                             | Required | Default
        agglom    | Number of cells to agglomerate          | no  | 0
    \endtable

SourceFiles
    randomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_randomDecomp_H
#define Foam_randomDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class randomDecomp Declaration
\*---------------------------------------------------------------------------*/

class randomDecomp
:
    public decompositionMethod
{
    // Private Data

        //- Number of cells to agglomerate per random value (default: 0)
        label agglom_;


    // Private Member Functions

        //- Random distribution on the 0-nCells interval
        labelList randomMap(const label nCells) const;

public:

    // Generated Methods

        //- No copy construct
        randomDecomp(const randomDecomp&) = delete;

        //- No copy assignment
        void operator=(const randomDecomp&) = delete;


    //- Runtime type information
    TypeNameNoDebug("random");


    // Constructors

        //- Construct with number of domains (no coefficients or constraints)
        explicit randomDecomp(const label numDomains);

        //- Construct for decomposition dictionary and optional region name
        explicit randomDecomp
        (
            const dictionary& decompDict,
            const word& regionName = "",
            int select = selectionType::DEFAULT
        );


    //- Destructor
    virtual ~randomDecomp() = default;


    // Member Functions

        //- Does not care about proc boundaries
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights_unused = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc_unused = pointField::null(),
            const scalarField& cWeights_unused = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,
            const pointField& cc_unused = pointField::null(),
            const scalarField& cWeights_unused = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc_unused = pointField::null(),
            const scalarField& cWeights_unused = scalarField::null()
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

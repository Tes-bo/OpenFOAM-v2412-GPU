/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::manualDecomp

Description
    Decompose based on cell-to-processor association in a file

    Method coefficients:
    \table
        Property  | Description                           | Required | Default
        dataFile  | filename of cell decomposition data   | yes      |
    \endtable

SourceFiles
    manualDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_manualDecomp_H
#define Foam_manualDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class manualDecomp Declaration
\*---------------------------------------------------------------------------*/

class manualDecomp
:
    public decompositionMethod
{
    // Private Data

        fileName dataFile_;

public:

    // Generated Methods

        //- No copy construct
        manualDecomp(const manualDecomp&) = delete;

        //- No copy assignment
        void operator=(const manualDecomp&) = delete;


    //- Runtime type information
    TypeName("manual");


    // Constructors

        //- Construct given decomposition dictionary and optional region name
        explicit manualDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~manualDecomp() = default;


    // Member Functions

        //- Manual decompose does not care about proc boundaries.
        //- Is all up to the user.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number.
        //  Use the mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
        ) const
        {
            NotImplemented;
            return labelList();
        }

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights = scalarField::null()
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

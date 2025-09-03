/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::structuredDecomp

Description
    Walk out decomposition of patch cells mesh - selectable as \c structured.

SourceFiles
    structuredDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_structuredDecomp_H
#define Foam_structuredDecomp_H

#include "decompositionMethod.H"
#include "wordRes.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class structuredDecomp Declaration
\*---------------------------------------------------------------------------*/

class structuredDecomp
:
    public decompositionMethod
{
    // Private Data

        dictionary methodDict_;

        wordRes patches_;

        autoPtr<decompositionMethod> method_;


public:

    // Generated Methods

        //- No copy construct
        structuredDecomp(const structuredDecomp&) = delete;

        //- No copy assignment
        void operator=(const structuredDecomp&) = delete;


    //- Runtime type information
    TypeName("structured");


    // Constructors

        //- Construct given decomposition dictionary. Region ignored
        explicit structuredDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~structuredDecomp() = default;


    // Member Functions

        //- Is method parallel aware
        //  (i.e. does it synchronize domains across proc boundaries)
        virtual bool parallelAware() const;

        //- Return for every coordinate the wanted processor number.
        //  Use the mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& points,
            const scalarField& pointWeights = scalarField::null()
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

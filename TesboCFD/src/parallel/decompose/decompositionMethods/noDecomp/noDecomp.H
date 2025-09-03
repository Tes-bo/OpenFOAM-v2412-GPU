/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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
    Foam::noDecomp

Description
    A dummy decomposition method, selected as \c none.

    Method coefficients: \a none

SourceFiles
    noDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_noDecomp_H
#define Foam_noDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class noDecomp Declaration
\*---------------------------------------------------------------------------*/

class noDecomp
:
    public decompositionMethod
{
public:

    // Generated Methods

        //- No copy construct
        noDecomp(const noDecomp&) = delete;

        //- No copy assignment
        void operator=(const noDecomp&) = delete;


    //- Runtime type information
    TypeNameNoDebug("none");


    // Constructors

        //- Construct with number of domains (ignored)
        explicit noDecomp(const label numDomains);

        //- Construct given decomposition dictionary and optional region name
        explicit noDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~noDecomp() = default;


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

        //- Return for every cell the current local processor rank.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc_unused = pointField::null(),
            const scalarField& cWeights_unused = scalarField::null()
        ) const;

        //- Return for every cell the current local processor rank.
        virtual labelList decompose
        (
            const CompactListList<label>& globalCellCells,
            const pointField& cc_unused = pointField::null(),
            const scalarField& cWeights_unused = scalarField::null()
        ) const;

        //- Return for every cell the current local processor rank.
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

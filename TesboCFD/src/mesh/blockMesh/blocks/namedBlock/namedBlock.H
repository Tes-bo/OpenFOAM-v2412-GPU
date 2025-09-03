/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::blocks::namedBlock

Description
    Gives name to a block.

SourceFiles
    namedBlock.C

\*---------------------------------------------------------------------------*/

#ifndef blocks_namedBlock_H
#define blocks_namedBlock_H

#include "block.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blocks
{

/*---------------------------------------------------------------------------*\
                           Class namedBlock Declaration
\*---------------------------------------------------------------------------*/

class namedBlock
:
    public word,
    public block
{
public:

    //- Runtime type information
    TypeNameNoDebug("name");


    // Constructors

        //- Construct from Istream setting points list
        namedBlock
        (
            const dictionary& dict,
            const label index,
            const pointField& vertices,
            const blockEdgeList& edges,
            const blockFaceList& faces,
            Istream& is
        );


    //- Destructor
    virtual ~namedBlock() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blocks
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

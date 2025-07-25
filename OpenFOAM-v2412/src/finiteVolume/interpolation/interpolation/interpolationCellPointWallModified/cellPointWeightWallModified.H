/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::cellPointWeightWallModified

Description
    Foam::cellPointWeightWallModified

SourceFiles
    cellPointWeightWallModified.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cellPointWeightWallModified_H
#define Foam_cellPointWeightWallModified_H

#include "cellPointWeight.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cellPointWeightWallModified Declaration
\*---------------------------------------------------------------------------*/

class cellPointWeightWallModified
:
    public cellPointWeight
{
public:

    // Constructors

        //- Construct from components
        cellPointWeightWallModified
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli,
            const label facei = -1
        );


    // Member Functions

        //- True if face is on a wall path
        static bool onWall(const polyMesh& mesh, const label facei);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

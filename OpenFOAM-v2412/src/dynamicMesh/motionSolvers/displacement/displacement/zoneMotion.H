/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::zoneMotion

Description
    Intermediate class for handling "zonified" motion.
    Can select motion for the entire mesh, or restrict to a subset
    of points based on a specified cellSet or cellZone(s).

Usage
    \verbatim
    {
        cellSet     <name>;

        // and/or
        cellZone    <name or regex>;
    }
    \endverbatim

Note
    The cellZone can be specified as a wordRe and will select by name
    or group. The special name "none" is treated as ignored.

SourceFiles
    zoneMotion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_zoneMotion_H
#define Foam_zoneMotion_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class zoneMotion Declaration
\*---------------------------------------------------------------------------*/

class zoneMotion
{
     // Private Data

        //- Points to move when cell zones or sets are supplied
        labelList pointIDs_;

        //- Flag to indicate whether all cells should move
        bool moveAllCells_;


public:

    // Constructors

        //- Construct from dictionary
        zoneMotion(const dictionary& dict, const polyMesh& mesh);


    //- Destructor
    virtual ~zoneMotion() = default;


    // Member Functions

        //- The point ids (for cell set/zone subset)
        const labelList& pointIDs() const noexcept
        {
            return pointIDs_;
        }

        //- Move all cells?
        bool moveAllCells() const noexcept
        {
            return moveAllCells_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

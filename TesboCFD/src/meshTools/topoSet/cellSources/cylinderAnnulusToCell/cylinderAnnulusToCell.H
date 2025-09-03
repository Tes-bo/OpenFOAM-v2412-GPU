/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::cylinderAnnulusToCell

Group
    grpCellSources

Description
    A \c topoSetCellSource to select all cells with centres
    inside a given bounding cylinder annulus.

    Identical behaviour and parameters as Foam::cylinderToCell
    except the outer radius is mandatory in interactive mode.

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource
    - Foam::cylinderToCell

SourceFiles
    cylinderAnnulusToCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_cylinderAnnulusToCell_H
#define Foam_cylinderAnnulusToCell_H

#include "cylinderToCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cylinderAnnulusToCell Declaration
\*---------------------------------------------------------------------------*/

class cylinderAnnulusToCell
:
    public cylinderToCell
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

public:

    //- Runtime type information
    TypeName("cylinderAnnulusToCell");


    // Constructors

        //- Construct from components
        cylinderAnnulusToCell
        (
            const polyMesh& mesh,
            const point& point1,
            const point& point2,
            const scalar radius,
            const scalar innerRadius = 0
        );

        //- Construct from dictionary
        cylinderAnnulusToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cylinderAnnulusToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cylinderAnnulusToCell() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

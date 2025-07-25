/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::cellCentreSet

Description
    A sampleSet based on cell centres

    For a dictionary specification:
    \table
        Property | Description                             | Required | Default
        type     | cellCentre                              | yes      |
        axis     | x, y, z, xyz, distance                  | no       | xyz
        bounds   | limit with bounding box                 | no       |
    \endtable

SourceFiles
    cellCentreSet.C

\*---------------------------------------------------------------------------*/

#ifndef cellCentreSet_H
#define cellCentreSet_H

#include "sampledSet.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellCentreSet Declaration
\*---------------------------------------------------------------------------*/

class cellCentreSet
:
    public sampledSet
{
    // Private data

        //- Optional bounding box to restrict the geometry
        const boundBox bounds_;


    // Private Member Functions

        //- Generate samples
        void genSamples();


public:

    //- Runtime type information
    TypeName("cellCentre");


    // Constructors

        //- Construct from components
        cellCentreSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const boundBox& bbox = boundBox::null()
        );

        //- Construct from dictionary
        cellCentreSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellCentreSet() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

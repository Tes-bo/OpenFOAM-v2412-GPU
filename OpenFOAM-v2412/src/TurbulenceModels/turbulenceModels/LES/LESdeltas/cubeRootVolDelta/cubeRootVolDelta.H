/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::cubeRootVolDelta

Description
    Simple cube-root of cell volume delta used in LES models.

SourceFiles
    cubeRootVolDelta.C

\*---------------------------------------------------------------------------*/

#ifndef cubeRootVolDelta_H
#define cubeRootVolDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class cubeRootVolDelta Declaration
\*---------------------------------------------------------------------------*/

class cubeRootVolDelta
:
    public LESdelta
{
    // Private data

        scalar deltaCoeff_;


    // Private Member Functions

        //- No copy construct
        cubeRootVolDelta(const cubeRootVolDelta&) = delete;

        //- No copy assignment
        void operator=(const cubeRootVolDelta&) = delete;


public:

    //- Runtime type information
    TypeName("cubeRootVol");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        cubeRootVolDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~cubeRootVolDelta() = default;


    // Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        // Correct values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::maxDeltaxyzCubeRootLESDelta

Description
    Maximum delta between maxDeltaxyz and cubeRootVolDelta.

SourceFiles
    maxDeltaxyzCubeRootLESDelta.C

\*---------------------------------------------------------------------------*/

#ifndef maxDeltaxyzCubeRootLESDelta_H
#define maxDeltaxyzCubeRootLESDelta_H

#include "LESdelta.H"
#include "maxDeltaxyz.H"
#include "cubeRootVolDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                   Class maxDeltaxyzCubeRootLESDelta Declaration
\*---------------------------------------------------------------------------*/

class maxDeltaxyzCubeRootLESDelta
:
    public LESdelta
{

    // Private Member Functions

        //- No copy construct
        maxDeltaxyzCubeRootLESDelta
        (
            const maxDeltaxyzCubeRootLESDelta&
        ) = delete;

        //- No copy assignment
        void operator=(const maxDeltaxyzCubeRootLESDelta&) = delete;

        // Calculate the delta values
        void calcDelta();

        //- maxDeltaxyz model
        maxDeltaxyz maxDeltaxyz_;

        //- cubeRootVolDelta model
        cubeRootVolDelta cubeRootVolDelta_;


public:

    //- Runtime type information
    TypeName("maxDeltaxyzCubeRoot");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        maxDeltaxyzCubeRootLESDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~maxDeltaxyzCubeRootLESDelta() = default;


    // Member Functions

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

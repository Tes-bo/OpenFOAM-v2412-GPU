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
    Foam::functionObjects::runTimePostProcessingDummy

Group
    grpGraphicsFunctionObjects

Description
    Dummy implementation of Foam::functionObjects::runTimePostProcessing
    to report when the real version is unavailable.

SourceFiles
    runTimePostProcessingDummy.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostProcessingDummy_H
#define functionObjects_runTimePostProcessingDummy_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                 Class runTimePostProcessingDummy Declaration
\*---------------------------------------------------------------------------*/

class runTimePostProcessingDummy
:
    public functionObject::unavailableFunctionObject
{
public:

    //- Declare type-name, virtual type (with debug switch)
    TypeName("runTimePostProcessing");

    // Constructors

        //- Construct from dictionary
        runTimePostProcessingDummy
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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
    Foam::functionObjects::writeFreeSurface

Description
    A function object to write the control points on the free surface

SourceFiles
    writeFreeSurface.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_writeFreeSurface_H
#define Foam_functionObjects_writeFreeSurface_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class writeFreeSurface Declaration
\*---------------------------------------------------------------------------*/

class writeFreeSurface
:
    public functionObjects::fvMeshFunctionObject
{
    // Private Member Functions

        //- Write data unconditionally
        void writeData();


public:

    //- Runtime type information
    TypeName("writeFreeSurface");


    // Constructors

        //- Construct from Time and dictionary
        writeFreeSurface
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


        //- No copy construct
        writeFreeSurface(const writeFreeSurface&) = delete;

        //- No copy assignment
        void operator=(const writeFreeSurface&) = delete;


    //- Destructor
    virtual ~writeFreeSurface() = default;


    // Member Functions

        //- Called at each ++ or += of the time-loop
        virtual bool execute();

        //- No-op
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

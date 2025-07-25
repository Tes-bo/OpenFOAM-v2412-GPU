/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019, 2022 PCOpt/NTUA
    Copyright (C) 2013-2019, 2022 FOSS GP
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
    Foam::incompressible::RASVariables::laminar

Description

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef incompressible_RASVariables_laminar_H
#define incompressible_RASVariables_laminar_H

#include "RASModelVariables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASVariables
{

/*---------------------------------------------------------------------------*\
                       Class laminar Declaration
\*---------------------------------------------------------------------------*/

class laminar
:
    public RASModelVariables
{

public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Construct from components
        laminar
        (
            const fvMesh& mesh,
            const solverControl& SolverControl
        );


    //- Destructor
    virtual ~laminar() = default;


    // Member Functions

        // Bools to identify which turbulent fields are present
        inline virtual bool hasTMVar1() const
        {
            return false;
        }
        inline virtual bool hasTMVar2() const
        {
            return false;
        }
        inline virtual bool hasNut() const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASVariables
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

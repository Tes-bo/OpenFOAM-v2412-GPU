/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::faceCells

Description
    Smooth ATC in cells next to a set of patches supplied by type

SourceFiles
    faceCells.C

\*---------------------------------------------------------------------------*/

#ifndef faceCells_H
#define faceCells_H

#include "zeroATCcells.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class adjointTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class faceCells
:
    public zeroATCcells
{

private:

    // Private Member Functions

        //- No copy construct
        faceCells(const faceCells&) = delete;

        //- No copy assignment
        void operator=(const faceCells&) = delete;


public:

    //- Runtime type information
    TypeName("faceCells");


    // Constructors

        //- Construct from components
        faceCells
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~faceCells() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

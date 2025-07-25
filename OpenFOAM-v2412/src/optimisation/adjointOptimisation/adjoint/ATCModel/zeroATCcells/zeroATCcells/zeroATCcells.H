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
    Foam::zeroATCcells

Description
    Base class for selecting cells on which to zero the ATC term

SourceFiles
    zeroATCcells.C

\*---------------------------------------------------------------------------*/

#ifndef zeroATCcells_H
#define zeroATCcells_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class adjointTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class zeroATCcells
{
private:

    // Private Member Functions

        //- No copy construct
        zeroATCcells(const zeroATCcells&) = delete;

        //- No copy assignment
        void operator=(const zeroATCcells&) = delete;


protected:

    // Protected data

        const fvMesh& mesh_;
        wordList zeroATCPatches_;
        labelList zeroATCZones_;
        labelList zeroATCcells_;


public:

    //- Runtime type information
    TypeName("zeroATCcells");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            zeroATCcells,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

    // Constructors

        //- Construct from components
        zeroATCcells
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<zeroATCcells> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~zeroATCcells() = default;


    // Member Functions

        //- Get the zeroATCcells
        inline const labelList& getZeroATCcells()
        {
            return zeroATCcells_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

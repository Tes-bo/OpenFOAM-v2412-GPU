/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::dummyAgglomeration

Description
    Agglomerate without combining cells. Used for testing.

SourceFiles
    dummyAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dummyAgglomeration_H
#define Foam_dummyAgglomeration_H

#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dummyAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class dummyAgglomeration
:
    public GAMGAgglomeration
{
    // Private Data

        //- Preset number of levels
        label nLevels_;


    // Private Member Functions

        //- No copy construct
        dummyAgglomeration(const dummyAgglomeration&) = delete;

        //- No copy assignment
        void operator=(const dummyAgglomeration&) = delete;


public:

    //- Runtime type information
    TypeName("dummy");


    // Constructors

        //- Construct given mesh and controls
        dummyAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Agglomerate from a starting level. Starting level is usually 0
        //- (initial mesh) but sometimes >0 (restarting after processor
        //- agglomeration)
        virtual void agglomerate
        (
            const label nCellsInCoarsestLevel,
            const label startLevel,
            const scalarField& startFaceWeights,
            const bool doProcessorAgglomerate = true
        )
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

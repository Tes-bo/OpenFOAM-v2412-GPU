/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::pairGAMGAgglomeration

Description
    Agglomerate using the pair algorithm.

SourceFiles
    pairGAMGAgglomeration.C
    pairGAMGAgglomerate.C

\*---------------------------------------------------------------------------*/

#ifndef pairGAMGAgglomeration_H
#define pairGAMGAgglomeration_H

#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pairGAMGAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class pairGAMGAgglomeration
:
    public GAMGAgglomeration
{
    // Private data

        //- Number of levels to merge, 1 = don't merge, 2 = merge pairs etc.
        label mergeLevels_;

        //- Direction of cell loop for the current level
        static bool forward_;


protected:

    // Protected Member Functions

        //- No copy construct
        pairGAMGAgglomeration(const pairGAMGAgglomeration&) = delete;

        //- No copy assignment
        void operator=(const pairGAMGAgglomeration&) = delete;


public:

    //- Runtime type information
    TypeName("pair");


    // Constructors

        //- Construct given mesh and controls
        pairGAMGAgglomeration
        (
            const lduMesh& mesh,
            const dictionary& controlDict
        );

        //- Calculate and return agglomeration
        static tmp<labelField> agglomerate
        (
            label& nCoarseCells,
            const lduAddressing& fineMatrixAddressing,
            const scalarField& faceWeights
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
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

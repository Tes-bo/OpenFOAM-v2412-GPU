/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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
    Foam::masterCoarsestGAMGProcAgglomeration

Description
    Processor agglomeration of GAMGAgglomerations.

    - by default agglomerates onto the master processor
    - optionally have multiple masters through the nProcessorsPerMaster
      parameter. e.g.

        p
        {
            solver                  GAMG;
            tolerance               1e-06;
            relTol                  0.1;
            smoother                GaussSeidel;

            nCellsInCoarsestLevel   10;

            processorAgglomerator   masterCoarsest;
            // Groups of 32 cores get combined so each 'master' gets
            // (roughly) 32*nCellsInCoarsestLevel cells.
            nProcessorsPerMaster    32;

            // Alternative : specify number of (equi-distributed) masters
            nMasters                2;

            // Restart local agglomeration after processor agglomeration. Used
            // as nCellsInCoarsestLevel. Only applicable for pair-wise local
            // agglomeration methods.
            nCellsInMasterLevel     1;
        }

SourceFiles
    masterCoarsestGAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_masterCoarsestGAMGProcAgglomeration_H
#define Foam_masterCoarsestGAMGProcAgglomeration_H

#include "GAMGProcAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class GAMGAgglomeration;

/*---------------------------------------------------------------------------*\
              Class masterCoarsestGAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class masterCoarsestGAMGProcAgglomeration
:
    public GAMGProcAgglomeration
{
    // Private Data

        label nProcessorsPerMaster_;

        const label nCellsInMasterLevel_;

        DynamicList<label> comms_;


    // Private Member Functions

        //- No copy construct
        masterCoarsestGAMGProcAgglomeration
        (
            const masterCoarsestGAMGProcAgglomeration&
        ) = delete;

        //- No copy assignment
        void operator=(const masterCoarsestGAMGProcAgglomeration&) = delete;


public:

    //- Runtime type information
    TypeName("masterCoarsest");


    // Constructors

        //- Construct given agglomerator and controls
        masterCoarsestGAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor
    virtual ~masterCoarsestGAMGProcAgglomeration();


    // Member Functions

        //- Modify agglomeration. Return true if modified
        virtual bool agglomerate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

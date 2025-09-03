/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::GAMGProcAgglomeration

Description
    Processor agglomeration of GAMGAgglomerations.

SourceFiles
    GAMGProcAgglomeration.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_GAMGProcAgglomeration_H
#define Foam_GAMGProcAgglomeration_H

#include "runTimeSelectionTables.H"
#include "labelList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class GAMGAgglomeration;
class lduMesh;
class lduPrimitiveMesh;

/*---------------------------------------------------------------------------*\
                    Class GAMGProcAgglomeration Declaration
\*---------------------------------------------------------------------------*/

class GAMGProcAgglomeration
{
protected:

    // Protected Data

        //- Reference to agglomeration
        GAMGAgglomeration& agglom_;

        //- Allocated communicators
        DynamicList<label> comms_;


    // Protected Member Functions

        //- Debug: write agglomeration info
        void printStats(Ostream& os, GAMGAgglomeration& agglom) const;

        //- Agglomerate a level. Return true if anything has changed
        bool agglomerate
        (
            const label fineLevelIndex,
            const labelList& procAgglomMap,
            const labelList& masterProcs,
            const List<label>& agglomProcIDs,
            const label procAgglomComm
        );

        //- Debug: calculate global cell-cells
        static labelListList globalCellCells(const lduMesh&);

        //- Clear/free allocated communicators
        void clearCommunicators();


public:

    //- Runtime type information
    TypeName("GAMGProcAgglomeration");


    // Declare run-time constructor selection tables

        //- Runtime selection table for pure geometric agglomerators
        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGProcAgglomeration,
            GAMGAgglomeration,
            (
                GAMGAgglomeration& agglom,
                const dictionary& controlDict
            ),
            (
                agglom,
                controlDict
            )
        );


    //- No copy construct
    GAMGProcAgglomeration(const GAMGProcAgglomeration&) = delete;

    //- No copy assignment
    void operator=(const GAMGProcAgglomeration&) = delete;


    // Constructors

        //- Construct given agglomerator and controls
        GAMGProcAgglomeration
        (
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    // Selectors

        //- Return the selected agglomerator
        static autoPtr<GAMGProcAgglomeration> New
        (
            const word& type,
            GAMGAgglomeration& agglom,
            const dictionary& controlDict
        );


    //- Destructor. Clears allocated communicators
    virtual ~GAMGProcAgglomeration();


    // Member Functions

        //- Modify agglomeration. \returns True if modified
        virtual bool agglomerate() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

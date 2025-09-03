/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::fileOperations::hostUncollatedFileOperation

Description
    Version of masterUncollated with host-based IO ranks

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileOperations_hostUncollatedFileOperation_H
#define Foam_fileOperations_hostUncollatedFileOperation_H

#include "masterUncollatedFileOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                 Class hostUncollatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class hostUncollatedFileOperation
:
    public masterUncollatedFileOperation
{
    // Private Data

        //- Communicator allocated/managed by us
        mutable label managedComm_;


    // Private Member Functions

        //- Any initialisation steps after constructing
        void init(bool verbose);

public:

        //- Runtime type information
        TypeName("hostUncollated");


    // Constructors

        //- Default construct
        explicit hostUncollatedFileOperation(bool verbose = false);

        //- Construct from communicator with specified io-ranks
        explicit hostUncollatedFileOperation
        (
            const Tuple2<label, labelList>& commAndIORanks,
            const bool distributedRoots,
            bool verbose = false
        );


    //- Destructor
    virtual ~hostUncollatedFileOperation();


    // Member Functions

        //- Transfer ownership of communicator to this fileOperation.
        //- Use with caution
        virtual void storeComm() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

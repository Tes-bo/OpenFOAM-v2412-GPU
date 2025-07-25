/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::fileOperations::collatedFileOperation

Description
    Version of masterUncollatedFileOperation that collates regIOobjects
    into a container in the processors/ subdirectory.

    Uses threading if maxThreadFileBufferSize != 0.
        > 0 : Can use mpi inside thread to collect data if buffer is not
              large enough. Does need full thread support inside MPI.

        < 0 : special : -maxThreadFileBufferSize is guaranteed large enough
              for all writing. Initialises MPI without thread support.

See also
    masterUncollatedFileOperation

SourceFiles
    collatedFileOperation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileOperations_collatedFileOperation_H
#define Foam_fileOperations_collatedFileOperation_H

#include "masterUncollatedFileOperation.H"
#include "fileOperationInitialise.H"
#include "OFstreamCollator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                    Class collatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class collatedFileOperation
:
    public masterUncollatedFileOperation
{
    // Private Data

        //- Communicator allocated/managed by us
        mutable label managedComm_;


    // Private Member Functions

        //- Any initialisation steps after constructing
        void init(bool verbose);


protected:

    // Protected Data

        //- Threaded writer
        mutable OFstreamCollator writer_;


    // Protected Member Functions

        //- Print banner information, optionally with io ranks
        void printBanner(const bool withRanks = false) const;

        //- Append to processorsNN/ file
        bool appendObject
        (
            const regIOobject& io,
            const fileName& pathName,
            IOstreamOption streamOpt
        ) const;


public:

    //- Runtime type information
    TypeName("collated");


    // Static Data

        //- Max size of thread buffer size. This is the overall size of
        //  all files to be written. Starts blocking if not enough size.
        //  Read as float to enable easy specification of large sizes.
        static float maxThreadFileBufferSize;


    // Constructors

        //- Default construct
        explicit collatedFileOperation(bool verbose = false);

        //- Construct from communicator with specified io-ranks
        explicit collatedFileOperation
        (
            const Tuple2<label, labelList>& commAndIORanks,
            const bool distributedRoots,
            bool verbose = false
        );


    //- Destructor
    virtual ~collatedFileOperation();


    // Member Functions

        //- Transfer ownership of communicator to this fileOperation.
        //- Use with caution
        virtual void storeComm() const;


        // (reg)IOobject functionality

            //- Generate disk file name for object. Opposite of filePath.
            virtual fileName objectPath
            (
                const IOobject& io,
                const word& typeName
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state.
            virtual bool writeObject
            (
                const regIOobject&,
                IOstreamOption streamOpt = IOstreamOption(),
                const bool writeOnProc = true
            ) const;


        // Other

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Actual name of processors dir
            virtual word processorsDir(const IOobject&) const;

            //- Actual name of processors dir
            virtual word processorsDir(const fileName&) const;
};


/*---------------------------------------------------------------------------*\
               Class fileOperationInitialise_collated Declaration
\*---------------------------------------------------------------------------*/

//- A fileOperation initialiser for collated file handlers.
//- Requires threading for non-zero maxThreadFileBufferSize.
class fileOperationInitialise_collated
:
    public fileOperationInitialise
{
public:

    // Constructors

        //- Construct from components
        fileOperationInitialise_collated(int& argc, char**& argv)
        :
            fileOperationInitialise(argc, argv)
        {}


    //- Destructor
    virtual ~fileOperationInitialise_collated() = default;


    // Member Functions

        //- The (MPI) threading requirement depends on buffering
        virtual bool needsThreading() const
        {
            return (collatedFileOperation::maxThreadFileBufferSize > 0);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

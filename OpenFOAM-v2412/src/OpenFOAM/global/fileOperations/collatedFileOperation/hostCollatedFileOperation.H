/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::fileOperations::hostCollatedFileOperation

Description
    Version of collatedFileOperation with multiple read/write ranks.

    In parallel it will assume ranks are sorted according to hostname
    and the lowest rank per hostname will be the IO rank. The output directories
    will get a unique name processors\<N\>_\<low\>-\<high\> where N is the
    overall number of processors and low and high is the range of ranks
    contained in the files. Each of these subsets uses its own communicator.

    Instead of using the hostnames the IO ranks can be assigned using the
    FOAM_IORANKS environment variable (also when running non-parallel), e.g.
    when decomposing into 4:

    \verbatim
        FOAM_IORANKS='(0 2)' decomposePar -fileHandler hostCollated
        FOAM_IORANKS='0,2'   decomposePar -fileHandler hostCollated
        FOAM_IORANKS='0 2'   decomposePar -fileHandler hostCollated
    \endverbatim

    will generate

    \verbatim
        processors4_0-1/  : containing data for processors 0 to 1
        processors4_2-3/  : containing data for processors 2 to 3
    \endverbatim

Environment
    - \c FOAM_ENV : list of io-ranks as plain space or comma separated
      list or as an OpenFOAM formatted list. Eg, '(0 4 8)'

See also
    Foam::collatedFileOperation

SourceFiles
    hostCollatedFileOperation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileOperations_hostCollatedFileOperation_H
#define Foam_fileOperations_hostCollatedFileOperation_H

#include "collatedFileOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                 Class hostCollatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class hostCollatedFileOperation
:
    public collatedFileOperation
{
    // Private Data

        //- Communicator allocated/managed by us
        mutable label managedComm_;


    // Private Member Functions

        //- Any initialisation steps after constructing
        void init(bool verbose);


public:

        //- Runtime type information
        TypeName("hostCollated");


    // Constructors

        //- Default construct
        explicit hostCollatedFileOperation(bool verbose = false);

        //- Construct from communicator with specified io-ranks
        explicit hostCollatedFileOperation
        (
            const Tuple2<label, labelList>& commAndIORanks,
            const bool distributedRoots,
            bool verbose = false
        );


    //- Destructor
    virtual ~hostCollatedFileOperation();


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

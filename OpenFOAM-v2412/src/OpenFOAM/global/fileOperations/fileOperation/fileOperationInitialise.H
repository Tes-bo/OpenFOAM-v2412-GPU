/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::fileOperationInitialise

Description
    General fileOperation initialiser.
    Handles \c -ioRanks option, using it to set the FOAM_IORANKS environment
    variable.

\*---------------------------------------------------------------------------*/

#ifndef Foam_fileOperations_fileOperationInitialise_H
#define Foam_fileOperations_fileOperationInitialise_H

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                   Class fileOperationInitialise Declaration
\*---------------------------------------------------------------------------*/

class fileOperationInitialise
{
public:

        //- Runtime type information
        TypeName("fileOperationInitialise");


    // Constructors

        //- Construct from components
        fileOperationInitialise(int& argc, char**& argv);


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fileOperationInitialise,
            word,
            (
                int& argc, char**& argv
            ),
            (argc, argv)
        );


    // Selectors

        //- Select initialisation type
        static autoPtr<fileOperationInitialise> New
        (
            const word& type, int& argc, char**& argv
        );


    //- Destructor
    virtual ~fileOperationInitialise() = default;


    // Member Functions

        //- Requires (MPI) threading?
        virtual bool needsThreading() const = 0;
};


/*---------------------------------------------------------------------------*\
              Class fileOperationInitialise_unthreaded Declaration
\*---------------------------------------------------------------------------*/

//- A fileOperation initialiser for unthreaded file handlers.
class fileOperationInitialise_unthreaded
:
    public fileOperationInitialise
{
public:

    // Constructors

        //- Construct from components
        fileOperationInitialise_unthreaded(int& argc, char**& argv)
        :
            fileOperationInitialise(argc, argv)
        {}


    //- Destructor
    virtual ~fileOperationInitialise_unthreaded() = default;


    // Member Functions

        //- No (MPI) threading required
        virtual bool needsThreading() const { return false; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

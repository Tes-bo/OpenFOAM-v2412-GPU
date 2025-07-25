/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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
    Foam::singleRun

Description
    Update design variables using steepest descent

SourceFiles
    singleRun.C

\*---------------------------------------------------------------------------*/

#ifndef singleRun_H
#define singleRun_H

#include "optimisationManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class singleRun Declaration
\*---------------------------------------------------------------------------*/

class singleRun
:
    public optimisationManager
{
protected:

    // Protected data
    label cycles_;


private:

    // Private Member Functions

        //- No copy construct
        singleRun(const singleRun&) = delete;

        //- No copy assignment
        void operator=(const singleRun&) = delete;


public:

    //- Runtime type information
    TypeName("singleRun");


    // Constructors

        //- Construct from components
        singleRun(fvMesh& mesh);


    //- Destructor
    virtual ~singleRun() = default;


    // Member Functions

        //- Prefix increment
        virtual optimisationManager& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual optimisationManager& operator++(int);

        //- Return true if end of optimisation run.
        //  Also, updates the design variables if needed
        virtual bool checkEndOfLoopAndUpdate();

        //- Return true if end of optimisation run
        virtual bool end();

        //- Whether to update the design variables
        virtual bool update();

        //- Update design variables.
        //  Does nothing
        virtual void updateDesignVariables();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

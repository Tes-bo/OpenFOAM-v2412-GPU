/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::steadyOptimisation

Description
    Iterate the optimisation cycles. For steady state opt, this coinsides
    with evolving Time

SourceFiles
    steadyOptimisation.C

\*---------------------------------------------------------------------------*/

#ifndef steadyOptimisation_H
#define steadyOptimisation_H

#include "optimisationManager.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class steadyOptimisation Declaration
\*---------------------------------------------------------------------------*/

class steadyOptimisation
:
    public optimisationManager
{
private:

    // Private Member Functions

        //- No copy construct
        steadyOptimisation(const steadyOptimisation&) = delete;

        //- No copy assignment
        void operator=(const steadyOptimisation&) = delete;


protected:

    // Protected Member Functions

        //- Update design variables using a line-search
        void lineSearchUpdate();

        //- Update design variables using a fixed step
        void fixedStepUpdate();


public:

    //- Runtime type information
    TypeName("steadyOptimisation");


    // Constructors

        //- Construct from components
        steadyOptimisation(fvMesh& mesh);


    // Destructor
    virtual ~steadyOptimisation() = default;


    // Member Functions

        //- Prefix increment
        virtual optimisationManager& operator++();

        //- Postfix increment, this is identical to the prefix increment
        virtual optimisationManager& operator++(int);

        //- Return true if end of optimisation run
        //  Also, updates the design variables if needed
        virtual bool checkEndOfLoopAndUpdate();

        //- Return true if end of optimisation run
        virtual bool end();

        //- Whether to update the design variables
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

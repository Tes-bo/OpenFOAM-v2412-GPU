/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::simpleControl

Description
    SIMPLE control class to supply convergence information/checks for
    the SIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef simpleControl_H
#define simpleControl_H

#include "solutionControl.H"

//- Declare that pimpleControl will be used
#define SIMPLE_CONTROL

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class simpleControl Declaration
\*---------------------------------------------------------------------------*/

class simpleControl
:
    public solutionControl
{

protected:

    // Protected Data

        //- Initialised flag
        bool initialised_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        bool read();

        //- Return true if all convergence checks are satisfied
        bool criteriaSatisfied();


private:

    // Private member functions

        //- No copy construct
        simpleControl(const simpleControl&) = delete;

        //- No copy assignment
        void operator=(const simpleControl&) = delete;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("simpleControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        simpleControl
        (
            fvMesh& mesh,
            const word& dictName = "SIMPLE",
            const bool verbose = true
        );


    //- Destructor
    virtual ~simpleControl() = default;


    // Member Functions

        // Solution control

            //- SIMPLE loop
            virtual bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

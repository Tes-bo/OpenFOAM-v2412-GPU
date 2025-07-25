/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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
    Foam::dynamicTopODesignVariables

Description
    Design variables for topology optimisation problems.
    Active design variables are added gradually, to help the optimiser from
    being trapped to a sub-optimal solution track in the first optimisation
    cycles.

SourceFiles
    dynamicTopODesignVariables.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicTopODesignVariables_H
#define dynamicTopODesignVariables_H

#include "topODesignVariables.H"
#include "marchingCells.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class dynamicTopODesignVariables Declaration
\*---------------------------------------------------------------------------*/

class dynamicTopODesignVariables
:
    public topODesignVariables
{
protected:

    // Protected data

        //- Mechanism for gradually activating design variables
        marchingCells marchCells_;

        //- Number of times the design variables have evolved
        label evolvedCount_;


    // Protected Member Functions

        //- Set active design variables
        virtual void setActiveDesignVariables
        (
            const label fluidID = 0,
            const bool activeIO = false
        );


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicTopODesignVariables
        (
            const dynamicTopODesignVariables&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const dynamicTopODesignVariables&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicTopO");


    // Constructors

        //- Construct from dictionary
        dynamicTopODesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct from dictionary and size
        dynamicTopODesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict,
            const label size
        );


    // Selectors

        //- Construct and return the selected design variables
        static autoPtr<dynamicTopODesignVariables> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~dynamicTopODesignVariables() = default;


    // Member Functions

        //- Update the active design variables
        virtual void evolveNumber();

        //- The writeData function required by the regIOobject write operation
        virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

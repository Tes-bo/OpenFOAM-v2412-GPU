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
    Foam::ArmijoConditions

Description
    Class satisfying the Armijo line search conditions

SourceFiles
    ArmijoConditions.C

\*---------------------------------------------------------------------------*/

#ifndef ArmijoConditions_H
#define ArmijoConditions_H

#include "lineSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ArmijoConditions Declaration
\*---------------------------------------------------------------------------*/

class ArmijoConditions
:
    public lineSearch
{
protected:

    // Protected data

        //- Multiplier of the merit function reduction computed using
        //- a first-order Taylor expansion
        scalar c1_;


private:

    // Private Member Functions

        //- No copy construct
        ArmijoConditions(const ArmijoConditions&) = delete;

        //- No copy assignment
        void operator=(const ArmijoConditions&) = delete;


public:

    //- Runtime type information
    TypeName("ArmijoConditions");


    // Constructors

        //- Construct from components
        ArmijoConditions
        (
            const dictionary& dict,
            const Time& time,
            updateMethod& UpdatheMethod
        );


    // Destructor
    virtual ~ArmijoConditions() = default;


    // Member Functions

       //Return the correction of the design variables
       virtual bool converged();

       //Update step in given direction
       virtual void updateStep();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::constraintProjection

Description
    Update design variables using Rosen's projection method.

    Can handle multiple constraints.
    Optionally includes a non-linearity correction

SourceFiles
    constraintProjection.C

\*---------------------------------------------------------------------------*/

#ifndef constraintProjection_H
#define constraintProjection_H

#include "constrainedOptimisationMethod.H"
#include "updateMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class constraintProjection Declaration
\*---------------------------------------------------------------------------*/

class constraintProjection
:
    public constrainedOptimisationMethod,
    public updateMethod
{
protected:

    // Protected data

        //- Correct for non-linearities
        bool useCorrection_;

        //- Weight of the aggregated constraint values in the merit function
        scalar delta_;


private:

    // Private Member Functions

        //- No copy construct
        constraintProjection(const constraintProjection&) = delete;

        //- No copy assignment
        void operator=(const constraintProjection&) = delete;


public:

    //- Runtime type information
    TypeName("constraintProjection");


    // Constructors

        //- Construct from components
        constraintProjection
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~constraintProjection() = default;


    // Member Functions

        //- Compute design variables correction
        void computeCorrection();

        //- Compute merit function. Could be different than the objective
        //- in the presence of constraints
        virtual scalar computeMeritFunction();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

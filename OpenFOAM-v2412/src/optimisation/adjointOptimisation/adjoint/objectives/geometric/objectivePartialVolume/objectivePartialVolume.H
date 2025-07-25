/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::objectives::objectivePartialVolume

Description

SourceFiles
    objectivePartialVolume.C

\*---------------------------------------------------------------------------*/

#ifndef objectivePartialVolume_H
#define objectivePartialVolume_H

#include "objectiveGeometric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                   Class objectivePartialVolume Declaration
\*---------------------------------------------------------------------------*/

class objectivePartialVolume
:
    public objectiveGeometric
{
    // Private data

        scalar initVol_;
        labelHashSet objectivePatches_;


public:

    //- Runtime type information
    TypeName("partialVolume");


    // Constructors

        //- from components
        objectivePartialVolume
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectivePartialVolume() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update d (x) / db multiplier. Surface and volume-based sensitivity
        //- term
        virtual void update_dxdbDirectMultiplier();

        //- Update d (normal dS) / db multiplier. Surface and volume-based
        //- sensitivity term
        virtual void update_dSdbMultiplier();

        //- Write initial volume for continuation
        virtual bool writeData(Ostream& os) const;

        // Helper write functions

            //- Write headers for additional columns
            virtual void addHeaderInfo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

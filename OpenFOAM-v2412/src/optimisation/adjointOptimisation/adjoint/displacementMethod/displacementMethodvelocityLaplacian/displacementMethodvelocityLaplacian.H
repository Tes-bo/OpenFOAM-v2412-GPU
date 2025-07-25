/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::displacementMethodvelocityLaplacian

Description
    Wrapper class for the velocityLaplacian motion solver

SourceFiles
    displacementMethodvelocityLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef displacementMethodvelocityLaplacian_H
#define displacementMethodvelocityLaplacian_H

#include "displacementMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class displacementMethodvelocityLaplacian Declaration
\*---------------------------------------------------------------------------*/

class displacementMethodvelocityLaplacian
:
    public displacementMethod
{
protected:

    // Protected data

        pointVectorField& pointMotionU_;

        volVectorField& cellMotionU_;

        bool resetFields_;


private:

    // Private Member Functions

        //- No copy construct
        displacementMethodvelocityLaplacian
        (
            const displacementMethodvelocityLaplacian&
        ) = delete;

        //- No copy assignment
        void operator=(const displacementMethodvelocityLaplacian&) = delete;


public:

    //- Runtime type information
    TypeName("velocityLaplacian");


    // Constructors

        //- Construct from components
        displacementMethodvelocityLaplacian
        (
            fvMesh& mesh,
            const labelList& patchIDs
        );


    //- Destructor
    virtual ~displacementMethodvelocityLaplacian() = default;


    // Member Functions

       //- Set motion filed related to model based on given motion
       void setMotionField(const pointVectorField& pointMovement);

       //- Set motion filed related to model based on given motion
       void setMotionField(const volVectorField& cellMovement);

       //- Set control field as a vectorField. For methods working with
       //- parameters (RBF etc)
       void setControlField(const vectorField& controlField);

       //- Set control field as a vectorField. For methods working with
       //- parameters (RBF etc)
       void setControlField(const scalarField& controlField);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

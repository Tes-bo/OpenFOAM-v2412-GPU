/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
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
    Foam::controlPointsDefinition

Description
    Abstract base for selecting how to construct the control points
    of a volumetric B-Splines morpher

SourceFiles
    controlPointsDefinition.C

\*---------------------------------------------------------------------------*/

#ifndef controlPointsDefinition_H
#define controlPointsDefinition_H

#include "NURBS3DVolume.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class controlPointsDefinition Declaration
\*---------------------------------------------------------------------------*/

class controlPointsDefinition
{
protected:

    // Protected data

        //- The volumetric B-Splines box
        NURBS3DVolume& box_;

        //- The volumetric B-Splines control points
        vectorField& cps_;


    // Protected Member Functions

        //- Compute control points
        virtual void computeControlPoints() = 0;

        //- Transform control points using the equivalent position, rotation
        //- and scale values from the Paraview clip filter
        void transformControlPoints
        (
            const vector& geometryMin,
            const vector& geometryMax
        );


private:

    // Private Member Functions

        //- No copy construct
        controlPointsDefinition(const controlPointsDefinition&) = delete;

        //- No copy assignment
        void operator=(const controlPointsDefinition&) = delete;


public:

    //- Runtime type information
    TypeName("controlPointsDefinition");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            controlPointsDefinition,
            dictionary,
            (
                NURBS3DVolume& box
            ),
            (box)
        );


    // Constructors

        //- Construct from underlaying NURBS box
        controlPointsDefinition(NURBS3DVolume& box);


    // Selectors

        //- Return a reference to the selected controlPointsDefinition model
        static autoPtr<controlPointsDefinition> New
        (
            NURBS3DVolume& box
        );


    //- Destructor
    virtual ~controlPointsDefinition() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

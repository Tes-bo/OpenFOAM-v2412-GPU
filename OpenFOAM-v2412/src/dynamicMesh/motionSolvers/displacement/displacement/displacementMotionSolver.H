/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::displacementMotionSolver

Description
    Virtual base class for displacement motion solver

    The boundary displacement is set as a boundary condition
    on the pointDisplacement pointVectorField.

SourceFiles
    displacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_displacementMotionSolver_H
#define Foam_displacementMotionSolver_H

#include "points0MotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class displacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementMotionSolver
:
    public points0MotionSolver
{
protected:

    // Protected Data

        //- Point motion field
        mutable pointVectorField pointDisplacement_;


    // Protected Member Functions

        //- No copy construct
        displacementMotionSolver(const displacementMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const displacementMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("displacementMotionSolver");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            displacementMotionSolver,
            displacement,
            (
                const polyMesh& mesh,
                const IOdictionary& dict,
                const pointVectorField& pointDisplacement,
                const pointIOField& points0
            ),
            (mesh, dict, pointDisplacement, points0)
        );


    // Selectors

        //- Select constructed from polyMesh, dictionary and components.
        //  If dictionary was registered this will 'steal' that registration.
        static autoPtr<displacementMotionSolver> New
        (
            const word& solverTypeName,
            const polyMesh&,
            const IOdictionary&,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    // Constructors

        //- Construct from mesh and dictionary
        displacementMotionSolver
        (
            const polyMesh&,
            const IOdictionary&,
            const word& type
        );

        //- Construct from mesh and dictionary
        displacementMotionSolver
        (
            const polyMesh&,
            const IOdictionary&,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0,
            const word& type
        );


    //- Destructor
    virtual ~displacementMotionSolver() = default;


    // Member Functions

        //- Return reference to the point motion displacement field
        pointVectorField& pointDisplacement() noexcept
        {
            return pointDisplacement_;
        }

        //- Return const reference to the point motion displacement field
        const pointVectorField& pointDisplacement() const noexcept
        {
            return pointDisplacement_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

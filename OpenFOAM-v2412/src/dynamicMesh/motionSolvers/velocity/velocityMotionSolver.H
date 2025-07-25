/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::velocityMotionSolver

Description
    Virtual base class for velocity motion solver

    The boundary displacement is set as a boundary condition
    on the pointMotionU pointVectorField.

SourceFiles
    velocityMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityMotionSolver_H
#define velocityMotionSolver_H

#include "motionSolver.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class velocityMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class velocityMotionSolver
:
    public motionSolver
{
protected:

    // Protected Data

        //- Point motion field
        mutable pointVectorField pointMotionU_;


    // Protected Member Functions

        //- No copy construct
        velocityMotionSolver(const velocityMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const velocityMotionSolver&) = delete;

public:

    //- Runtime type information
    TypeName("velocityMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        velocityMotionSolver
        (
            const polyMesh&,
            const IOdictionary&,
            const word& type
        );


    //- Destructor
    virtual ~velocityMotionSolver() = default;


    // Member Functions

        //- Return reference to the point motion velocity field
        pointVectorField& pointMotionU()
        {
            return pointMotionU_;
        }

        //- Return const reference to the point motion velocity field
        const pointVectorField& pointMotionU() const
        {
            return pointMotionU_;
        }

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

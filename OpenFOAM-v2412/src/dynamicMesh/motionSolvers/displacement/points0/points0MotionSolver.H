/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::points0MotionSolver

Description
    Virtual base class for displacement motion solvers, where the point
    motion is relative to a set of fixed points (points0).

SourceFiles
    points0MotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_points0MotionSolver_H
#define Foam_points0MotionSolver_H

#include "motionSolver.H"
#include "pointFields.H"
#include "pointIOField.H"
#include "zoneMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class points0MotionSolver Declaration
\*---------------------------------------------------------------------------*/

class points0MotionSolver
:
    public motionSolver,
    public zoneMotion
{
protected:

    // Protected Data

        //- Starting points
        pointIOField points0_;


    // Protected Member Functions

        //- No copy construct
        points0MotionSolver(const points0MotionSolver&) = delete;

        //- No copy assignment
        void operator=(const points0MotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("points0MotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        points0MotionSolver
        (
            const polyMesh&,
            const IOdictionary&,
            const word& type
        );

        //- Construct from mesh and dictionary
        points0MotionSolver
        (
            const polyMesh&,
            const IOdictionary&,
            const pointIOField& points0,
            const word& type
        );


    //- Destructor
    virtual ~points0MotionSolver() = default;


    // Member Functions

        //- Return IO object for points0
        static IOobject points0IO(const polyMesh& mesh);

        //- Return reference to the reference ('0') pointField
        pointField& points0() noexcept
        {
            return points0_;
        }

        //- Return reference to the reference ('0') pointField
        const pointField& points0() const noexcept
        {
            return points0_;
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

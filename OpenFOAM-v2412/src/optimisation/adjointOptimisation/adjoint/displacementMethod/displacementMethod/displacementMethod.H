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
    Foam::displacementMethod

Description
    Abstract base class for displacement methods, which are a set or
    wrapper classes allowing to change the driving force of mesh motion solvers

SourceFiles
    displacementMethod.C

\*---------------------------------------------------------------------------*/

#ifndef displacementMethod_H
#define displacementMethod_H

#include "fvm.H"
#include "fvc.H"
#include "runTimeSelectionTables.H"
#include "addToRunTimeSelectionTable.H"

#include "motionSolver.H"
#include "volPointInterpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class displacementMethod Declaration
\*---------------------------------------------------------------------------*/

class displacementMethod
{
protected:

    // Protected data

        fvMesh& mesh_;

        //- IDs of the patches to be moved
        labelList patchIDs_;

        autoPtr<motionSolver> motionPtr_;

        scalar maxDisplacement_;

        //- Whether the motion solver prefers a point of a vol field as input
        bool preferPointField_;


private:

    // Private Member Functions

        //- No copy construct
        displacementMethod(const displacementMethod&) = delete;

        //- No copy assignment
        void operator=(const displacementMethod&) = delete;


public:

    //- Runtime type information
    TypeName("displacementMethod");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            displacementMethod,
            dictionary,
            (
                fvMesh& mesh,
                const labelList& patchIDs
            ),
            (mesh, patchIDs)
        );


    // Constructors

        //- Construct from components
        displacementMethod
        (
            fvMesh& mesh,
            const labelList& patchIDs
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<displacementMethod> New
        (
            fvMesh& mesh,
            const labelList& patchIDs
        );


    //- Destructor
    virtual ~displacementMethod() = default;


    // Member Functions

        //- Whether the motion solver prefers a point of a vol field as input
        virtual bool preferPointField() const;

        //- Set motion filed related to model based on given motion
        virtual void setMotionField(const pointVectorField& pointMovement) = 0;

        //- Set motion filed related to model based on given motion
        virtual void setMotionField(const volVectorField& cellMovement) = 0;

        //- Set control field as a vectorField. For methods working with
        //- parameters (RBF etc)
        virtual void setControlField(const vectorField& controlField) = 0;

        //- Set control field as a vectorField. For methods working with
        //- parameters (RBF etc)
        virtual void setControlField(const scalarField& controlField) = 0;

        //- Bound control field in certain directions etc. For methods working
        //- with parameters (RBF etc)
        //- does nothing by default
        virtual void boundControlField(vectorField& controlField);

        //- Get access to motionSolver
        autoPtr<motionSolver>& getMotionSolver();

        //- Get max displacement
        scalar getMaxDisplacement() const;

        //- Set parametertised patch IDs
        inline void setPatchIDs(const labelList& patchIDs)
        {
            patchIDs_ = patchIDs;
        }

        //- Update mesh
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

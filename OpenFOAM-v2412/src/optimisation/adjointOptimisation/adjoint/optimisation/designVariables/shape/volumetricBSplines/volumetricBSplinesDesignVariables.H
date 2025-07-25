/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::volumetricBSplinesDesignVariables

Description
    Volumetric B-Splines design variables for shape optimisation

SourceFiles
    volumetricBSplinesDesignVariables.C

\*---------------------------------------------------------------------------*/

#ifndef volumetricBSplinesDesignVariables_H
#define volumetricBSplinesDesignVariables_H

#include "shapeDesignVariables.H"
#include "volBSplinesBase.H"
#include "morphingBoxConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class volumetricBSplinesDesignVariables Declaration
\*---------------------------------------------------------------------------*/

class volumetricBSplinesDesignVariables
:
    public shapeDesignVariables,
    public localIOdictionary
{
protected:

    // Protected data

        //- Reference to underlaying volumetric B-Splines morpher
        volBSplinesBase& volBSplinesBase_;

        //- Should the control points be non-overlapping
        bool nonOverlappingCPs_;

        //- Should the bounds of the non-overlapping control points be updated
        //- in each optimisation cycle
        bool updateBounds_;

        //- Constraint imposed on the movement of the design variables.
        //  Can be used to e.g. impose a uniform movement of the control points
        //  in one direction, etc.
        autoPtr<morphingBoxConstraint> constraint_;


    // Protected Member Functions

        //- Size of the active control points
        virtual label sensSize() const;

        //- Components of the active control points
        virtual const labelList& activeSensitivities() const;

        //- Set IDs of active design variables
        //  Might be different than what volBSplinesBase_ returns, if
        //  constraint_ changes the number of design variables
        void setActiveDesignVariables();

        //- Set control points based on current design variables values
        void designVariablesToControlPoints();

        //- Set the design variables based on the current control points
        void controlPointsToDesignVariables();

        //- Set the design variables based on the given control points
        void controlPointsToDesignVariables(const vectorField& controlPoints);

        //- Read bounds for design variables, if present
        void readBounds
        (
            autoPtr<scalar> lowerBoundPtr = nullptr,
            autoPtr<scalar> upperBoundPtr = nullptr
        );

        //- Read one set of bounds (lower, upper)
        void readBounds
        (
            autoPtr<scalarField>& bounds,
            const word& boundsName,
            const label sign
        );

        //- Write current bounds to file
        void writeBounds(const scalarField& bounds, const word& name) const;

        //- Set the field driving the displacement method.
        //  Can be either the movement of the control points or the boundary
        //  displacement, depending on the method used to move the mesh
        void setDisplacement(const vectorField& cpMovement);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        volumetricBSplinesDesignVariables
        (
            const volumetricBSplinesDesignVariables&
        ) = delete;

        //- Disallow default bitwise assignment
        void operator=(const volumetricBSplinesDesignVariables&) = delete;


public:

    //- Runtime type information
    TypeName("volumetricBSplines");


    // Constructors

        //- Construct from components
        volumetricBSplinesDesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~volumetricBSplinesDesignVariables() = default;


    // Member Functions

        //- Are control points non-overlapping
        inline bool nonOverlappingCPs() const
        {
            return nonOverlappingCPs_;
        }

        //- Are bounds to be updated after each optmisation cycle
        inline bool updateBounds() const
        {
            return updateBounds_;
        }

        //- Constant access to the volBSplinesBase object
        inline const volBSplinesBase& getVolBSplinesBase() const
        {
            return volBSplinesBase_;
        }

        //- Non-constant access to the volBSplinesBase object
        inline volBSplinesBase& getVolBSplinesBase()
        {
            return volBSplinesBase_;
        }

        //- Transform correction of design variables to control points movement
        tmp<vectorField> controlPointMovement
        (
            const scalarField& correction
        );

        //- Update design variables based on a given correction
        virtual void update(scalarField& correction);

        //- Reset to starting point of line search
        virtual void resetDesignVariables();

        //- Compute eta if not set in the first step
        virtual scalar computeEta(scalarField& correction);

        //- Whether to use global sum when computing matrix-vector products
        //- in update methods
        virtual bool globalSum() const;

        //- Assemble the sensitivity derivatives, by also applying possible
        //- constraints
        virtual tmp<scalarField> assembleSensitivities
        (
            adjointSensitivity& adjointSens
        );

        //- For design variables with a dynamic character (i.e. changing
        //  number), perform the evolution.
        //  Hijacked here to evolve the bounds of the design variables
        //  after the end of each optimisation cycle
        virtual void evolveNumber();

        //- Write fields to support continuation
        virtual bool writeData(Ostream& os) const;

        // Fields related to sensitivity computations

            //- Get dxdb for all mesh points
            virtual tmp<vectorField> dxdbVol
            (
                const label varID
            ) const;

            //- Get dxdb for given design variable and patch
            virtual tmp<vectorField> dxdbFace
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dndb for given design variable and patch
            virtual tmp<vectorField> dndb
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dSdb for given design variable and patch
            virtual tmp<vectorField> dSdb
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dCdb for given design variable.
            //  Used for FI-based sensitivities
            virtual tmp<volVectorField> dCdb(const label varID) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

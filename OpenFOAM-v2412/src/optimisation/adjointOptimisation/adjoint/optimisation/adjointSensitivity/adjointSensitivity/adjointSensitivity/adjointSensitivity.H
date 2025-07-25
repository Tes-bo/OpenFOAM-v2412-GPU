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
    Foam::adjointSensitivity

Description
    Abstract base class for adjoint-based sensitivities

SourceFiles
    adjointSensitivity.C

\*---------------------------------------------------------------------------*/

#ifndef adjointSensitivityIncompressible_H
#define adjointSensitivityIncompressible_H

#include "boundaryFieldsFwd.H"
#include "adjointEikonalSolver.H"
#include "adjointMeshMovementSolver.H"
#include "sensitivity.H"
#include "volFieldsFwd.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class adjointSolver;

/*---------------------------------------------------------------------------*\
                     Class adjointSensitivity Declaration
\*---------------------------------------------------------------------------*/

class adjointSensitivity
:
    public sensitivity
{
protected:


    // Protected data

        //- Reference to the underlaying adjoint solver
        adjointSolver& adjointSolver_;

        //- The sensitivity derivative values
        scalarField derivatives_;

        //- Append this word to files related to the sensitivities
        word suffix_;

        //- Include distance variation in sensitivity computations
        bool includeDistance_;

        //- Adjoint eikonal equation solver
        autoPtr<adjointEikonalSolver> eikonalSolver_;

        //- Adjoint grid displacement solver
        autoPtr<adjointMeshMovementSolver> adjointMeshMovementSolver_;

        // Fields to accumulated through the adjoint solver

            // Shape optimisation

                //- Multiplier of grad(dx/b)
                autoPtr<volTensorField> gradDxDbMult_;

                //- Multiplier of div(dx/db)
                autoPtr<scalarField> divDxDbMult_;

                //- Multiplier of face dx/db
                //  The term that multiplies the adjoint-related part of the
                //  sensitivities in the (E)SI approach
                autoPtr<boundaryVectorField> dxdbMult_;

                //- Multiplier of dSf/db
                autoPtr<boundaryVectorField> dSfdbMult_;

                //- Multiplier of dnf/db
                autoPtr<boundaryVectorField> dnfdbMult_;

                //- Multiplier of dCf/db, found in the objective function
                autoPtr<boundaryVectorField> dxdbDirectMult_;

                //- Multiplier of dx/db computed at points,
                //- found in the objective function
                autoPtr<pointBoundaryVectorField> pointDxDbDirectMult_;

                //- Multiplier of dx/db, coming from boundary conditions that
                //- depend on the geometry, like rotatingWallVelocity
                autoPtr<boundaryVectorField> bcDxDbMult_;

                //- dx/db multiplier coming from fvOptions
                autoPtr<vectorField> optionsDxDbMult_;


private:

    // Private Member Functions

        //- No copy construct
        adjointSensitivity(const adjointSensitivity&) = delete;

        //- No copy assignment
        void operator=(const adjointSensitivity&) = delete;


public:

    //- Runtime type information
    TypeName("adjointSensitivity");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            adjointSensitivity,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                adjointSolver& adjointSolver
            ),
            (
                mesh,
                dict,
                adjointSolver
            )
        );


    // Constructors

        //- Construct from components
        adjointSensitivity
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<adjointSensitivity> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~adjointSensitivity() = default;


    // Member Functions

        //- Read dictionary if changed
        virtual bool readDict(const dictionary& dict);

        //- Const access to adjoint solver
        inline const adjointSolver& getAdjointSolver() const
        {
            return adjointSolver_;
        }

        //- Non-const access to adjoint solver
        inline adjointSolver& getAdjointSolver()
        {
            return adjointSolver_;
        }

        //- Should the adjoint eikonal PDE should be solved
        inline bool includeDistance() const
        {
            return includeDistance_;
        }

        //- Return the adjoint eikonal solver
        inline autoPtr<adjointEikonalSolver>& getAdjointEikonalSolver()
        {
            return eikonalSolver_;
        }

        //- Return the adjoint eikonal solver
        inline autoPtr<adjointMeshMovementSolver>&
        getAdjointMeshMovementSolver()
        {
            return adjointMeshMovementSolver_;
        }

        //- Set suffix
        inline void setSuffix(const word& suffix)
        {
            suffix_ = suffix;
        }

        //- Get suffix
        inline const word& getSuffix() const
        {
            return suffix_;
        }

        //- Should the parameterization compute the internalField of dxdb
        virtual bool computeDxDbInternalField() const;

        //- Accumulate sensitivity integrands
        //  Corresponds to the flow and adjoint part of the sensitivities
        virtual void accumulateIntegrand(const scalar dt) = 0;

        //- Assemble sensitivities
        //  Adds the geometric part of the sensitivities
        virtual void assembleSensitivities
        (
            autoPtr<designVariables>& designVars
        );

        //- Calculates and returns sensitivity fields.
        //  Used with optimisation libraries
        virtual const scalarField& calculateSensitivities
        (
            autoPtr<designVariables>& designVars
        );

        //- Returns the sensitivity fields
        //  Assumes it has already been updated/computed
        const scalarField& getSensitivities() const;

        //- Zero sensitivity fields and their constituents
        virtual void clearSensitivities();

        //- Write sensitivity fields.
        //  If valid, copies boundaryFields to volFields and writes them.
        //  Virtual to be reimplemented by control points-based methods
        //  (Bezier, RBF) which do not need to write fields
        virtual void write(const word& baseName = word::null);

        // Access functions to multipliers

            // Shape optimisation

                inline const autoPtr<volTensorField>& gradDxDbMult() const;
                inline autoPtr<volTensorField>& gradDxDbMult();
                inline const autoPtr<scalarField>& divDxDbMult() const;
                inline const autoPtr<boundaryVectorField>& dxdbMult() const;
                inline const autoPtr<boundaryVectorField>& dSfdbMult() const;
                inline const autoPtr<boundaryVectorField>& dnfdbMult() const;
                inline const autoPtr<boundaryVectorField>&
                    dxdbDirectMult() const;
                inline const autoPtr<pointBoundaryVectorField>&
                    pointDxDbDirectMult() const;
                inline const autoPtr<boundaryVectorField>& bcDxDbMult() const;
                inline const autoPtr<vectorField>& optionsDxDbMult() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "adjointSensitivityI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

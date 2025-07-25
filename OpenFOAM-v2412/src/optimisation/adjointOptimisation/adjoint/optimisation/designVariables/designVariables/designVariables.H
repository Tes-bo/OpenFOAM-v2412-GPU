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
    Foam::designVariables

Description
    Abstract base class for defining design variables.

SourceFiles
    designVariables.C

\*---------------------------------------------------------------------------*/

#ifndef designVariables_H
#define designVariables_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "volFields.H"
#include "dictionary.H"
#include "primalSolver.H"
#include "adjointSolverManager.H"
#include "adjointSensitivity.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class adjointSensitivity;

/*---------------------------------------------------------------------------*\
                       Class designVariables Declaration
\*---------------------------------------------------------------------------*/

class designVariables
:
    public scalarField
{
protected:

    // Protected data

        fvMesh& mesh_;
        dictionary dict_;

        //- Which of the design variables will be updated
        labelList activeDesignVariables_;

        //- Copy of old design variables. Useful when performing line-search
        autoPtr<scalarField> oldDesignVariables_;

        //- Maximum design variables' change in the first optimisation cycle
        //  Used when eta is not used in updateMethod
        autoPtr<scalar> maxInitChange_;

        //- Lower bounds of the design variables
        autoPtr<scalarField> lowerBounds_;

        //- Upper bounds of the design variables
        autoPtr<scalarField> upperBounds_;


    // Protected Member Functions

        //- Read bounds for design variables, if present
        void readBounds
        (
            autoPtr<scalar> lowerBoundPtr = nullptr,
            autoPtr<scalar> upperBoundPtr = nullptr
        );


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        designVariables(const designVariables&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const designVariables&) = delete;


public:

    //- Runtime type information
    TypeName("designVariables");


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            designVariables,
            designVariables,
            (
                fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from dictionary
        designVariables
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct from dictionary and size
        designVariables
        (
            fvMesh& mesh,
            const dictionary& dict,
            const label size
        );


    // Selectors

        //- Return a reference to the selected design variables
        static autoPtr<designVariables> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~designVariables() = default;


    // Member Functions

        //- Read dictionary if changed
        virtual bool readDict(const dictionary& dict);

        //- Get the design variables
        //  Defaults to *this.
        //  Virtual for potential overriding from derived classes
        virtual const scalarField& getVars() const;

        //- Get the design variables
        //  Defaults to *this.
        //  Virtual for potential overriding from derived classes
        virtual scalarField& getVars();

        //- Update design variables based on a given correction
        //  Translates the scalarField of corrections to a meaningful
        //  update of the design variables
        virtual void update(scalarField& correction) = 0;

        //- Store design variables, as the starting point for line search
        virtual void storeDesignVariables();

        //- Reset to the starting point of line search
        virtual void resetDesignVariables();

        //- Compute eta if not set in the first step
        virtual scalar computeEta(scalarField& correction) = 0;

        //- Whether to use global sum when computing matrix-vector products
        //- in update methods
        //  Depends on whether the design variables are common for all
        //  processors (e.g. volumetric B-Splines control points) or distributed
        //  across them (e.g. topology optimisation)
        virtual bool globalSum() const = 0;

        //- Return list of active design variables
        inline const labelList& activeDesignVariables() const;

        //- Check whether the max. initial change of the design variables has
        //- been set
        inline bool isMaxInitChangeSet() const;

        //- Get maxInitChange
        inline const autoPtr<scalar>& getMaxInitChange();

        //- Set maxInitChange
        inline void setMaxInitChange(const scalar maxInitChange);

        //- Trigger the recomputation of eta by updateMethod
        inline virtual bool resetEta() const;

        //- Get min bounds for the design variables
        inline const autoPtr<scalarField>& lowerBounds() const;

        //- Get max bounds for the design variables
        inline const autoPtr<scalarField>& upperBounds() const;

        //- Get min bounds for the design variables
        inline scalarField& lowerBoundsRef();

        //- Get max bounds for the design variables
        inline scalarField& upperBoundsRef();

        //- Post process sensitivities if needed
        virtual void postProcessSens
        (
            scalarField& objectiveSens,
            PtrList<scalarField>& constraintSens,
            const wordList& adjointSolversNames,
            bool isMaster
        );

        //- Assemble sensitivity derivatives, by combining the part related
        //- to the primal and adjoint solution with the part related to the
        //- design variables
        virtual tmp<scalarField> assembleSensitivities
        (
            adjointSensitivity& sens
        ) = 0;

        //- For design variables with a dynamic character (i.e. changing
        //- number), perform the evolution
        virtual void evolveNumber();

        //- Set initial values of the design variables
        //  For design variables sets that need to be initialised after the
        //  construction of the primal fields.
        //  Does nothing in base
        virtual void setInitialValues();

        //- Add fvOptions depending on the design variables
        virtual void addFvOptions
        (
            const PtrList<primalSolver>& primalSolver,
            const PtrList<adjointSolverManager>& adjointSolverManagers
        );

        //- Design variables might add constraints related to themselves
        //- (e.g. linear combinations of the design variables)
        //- Return the values and gradients of these constraints
        virtual tmp<scalarField> constraintValues();
        virtual PtrList<scalarField> constraintDerivatives();

        //- Write useful quantities to files
        virtual void writeDesignVars();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "designVariablesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

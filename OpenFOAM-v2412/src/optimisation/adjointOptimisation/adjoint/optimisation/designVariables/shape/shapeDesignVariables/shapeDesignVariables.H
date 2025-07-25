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
    Foam::shapeDesignVariables

Description
    Abstract base class for defining design variables for shape optimisation.

SourceFiles
    shapeDesignVariables.C

\*---------------------------------------------------------------------------*/

#ifndef shapeDesignVariables_H
#define shapeDesignVariables_H

#include "designVariables.H"
#include "displacementMethod.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class shapeDesignVariables Declaration
\*---------------------------------------------------------------------------*/

class shapeDesignVariables
:
    public designVariables
{
protected:

    // Protected data

        //- Patches to be moved by the design variables
        labelHashSet parametertisedPatches_;

        //- Mesh movement mechanism
        autoPtr<displacementMethod> displMethodPtr_;

        //- Store old points. Useful for line search
        autoPtr<pointField> pointsInit_;

        //- Write the mesh points irrespective of whether this is a write time
        bool writeEachMesh_;

        // Auxiliary fields keeping track of the various sensitiity components

            //- Flow related term.
            //  Term including grad(dxdb) in the volume integrals of the FI
            //  formulation
            scalarField dxdbVolSens_;

            //- Flow related term.
            //  Main term in the E-SI formulation. Is the surface intergral
            //  emerging after performing the Gauss-divergence theorem on the
            //  FI-based sensitivities
            scalarField dxdbSurfSens_;

            //- Term depending on delta(n dS)/delta b
            scalarField dSdbSens_;

            //- Term depending on delta(n)/delta b
            scalarField dndbSens_;

            //- Term depending on delta(x)/delta b for objectives that directly
            //- depend on x
            scalarField dxdbDirectSens_;

            //- Term depending on delta(V)/delta b
            scalarField dVdbSens_;

            //- Term depending on distance differentiation
            scalarField distanceSens_;

            //- Term depending on fvOptions
            scalarField optionsSens_;

            //- Term depending on the differenation of boundary conditions
            scalarField bcSens_;

        //- Name of the sensitivity derivatives folder
        fileName derivativesFolder_;


    // Protected Member Functions

        //- Size of the sensitivity derivatives.
        //  Might be different than this->size() in some cases
        virtual label sensSize() const;

        //- Active variables for which to compute sensitivities
        //  Might be different than this->activeDesignVariables_ in some cases
        virtual const labelList& activeSensitivities() const;

        //- Compute dxdb at the mesh cell centers by solving a Laplace PDE
        virtual tmp<volVectorField> solveMeshMovementEqn
        (
            const label patchI,
            const label varID
        ) const;

        //- Allocate the fields assosiated with the computation of sensitivities
        //  Not allocated in the constructor since the size of the design
        //  variables is usually not known there
        void allocateSensFields();

        //- Zero the fields assosiated with the computation of sensitivities
        void zeroSensFields();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        shapeDesignVariables(const shapeDesignVariables&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const shapeDesignVariables&) = delete;


public:

    //- Runtime type information
    TypeName("shape");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            shapeDesignVariables,
            dictionary,
            (
                fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        shapeDesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Construct and return the selected shapeDesignVariables
        static autoPtr<shapeDesignVariables> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~shapeDesignVariables() = default;


    // Member Functions

        //- Read dictionary if changed
        virtual bool readDict(const dictionary& dict);

        //- Update design variables based on a given correction.
        //  Translates the scalarField of corrections to a meaningful
        //  update of the design variables
        virtual void update(scalarField& correction) = 0;

        //- Store design variables, as the starting point for line search
        virtual void storeDesignVariables();

        //- Reset to starting point of line search
        virtual void resetDesignVariables();

        //- Compute eta if not set in the first step
        virtual scalar computeEta(scalarField& correction) = 0;

        //- Whether to use global sum when computing matrix-vector products
        //- in update methods
        //  Depends on whether the design variables are common for all
        //  processors (e.g. volumetric B-Splines control points) or distributed
        //  across them (e.g. topology optimisation)
        virtual bool globalSum() const = 0;

        // Functions related to mesh movement

            //- Move mesh based on displacementMethod
            virtual void moveMesh();

            //- Patches affected by the parameterisation
            inline const labelHashSet& getPatches() const
            {
                return parametertisedPatches_;
            }

            //- Return displacementMethod
            inline autoPtr<displacementMethod>& returnDisplacementMethod()
            {
                return displMethodPtr_;
            }


        // Fields related to sensitivity computations

            //- Add part of sensitivity derivatives related to geometry
            //- variations
            virtual tmp<scalarField> assembleSensitivities
            (
                adjointSensitivity& adjointSens
            );

            //- Write final sensitivity derivatives to files
            virtual void writeSensitivities
            (
                const scalarField& sens,
                const adjointSensitivity& adjointSens
            );

            //- Get dxdb for all mesh points
            virtual tmp<vectorField> dxdbVol
            (
                const label varID
            ) const;

            //- Get dxdb for a given design variable and patch
            virtual tmp<vectorField> dxdbFace
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dndb for a given design variable and patch
            virtual tmp<vectorField> dndb
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dSdb for a given design variable and patch
            virtual tmp<vectorField> dSdb
            (
                const label patchI,
                const label varID
            ) const;

            //- Get dCdb for a given design variable.
            //  Used for FI-based sensitivities
            virtual tmp<volVectorField> dCdb(const label varID) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

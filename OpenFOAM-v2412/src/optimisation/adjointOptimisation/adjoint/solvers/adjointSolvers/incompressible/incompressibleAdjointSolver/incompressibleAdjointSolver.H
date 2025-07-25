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
    Foam::incompressibleAdjointSolver

Description
    Base class for incompressibleAdjoint solvers

\*---------------------------------------------------------------------------*/

#ifndef incompressibleAdjointSolver_H
#define incompressibleAdjointSolver_H

#include "adjointSolver.H"
#include "incompressibleVars.H"
#include "incompressibleAdjointVars.H"
#include "ATCModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class incompressibleAdjointSolver Declaration
\*---------------------------------------------------------------------------*/

class incompressibleAdjointSolver
:
    public adjointSolver
{
private:

    // Private Member Functions

        //- No copy construct
        incompressibleAdjointSolver
        (
            const incompressibleAdjointSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const incompressibleAdjointSolver&) = delete;


protected:

    // Protected data

        //- Primal variable set
        incompressibleVars& primalVars_;

        //- Adjoint Transpose Convection options
        autoPtr<ATCModel> ATCModel_;

        //- Auxiliary bool to avoid a potentially expensive
        //- part of the sensitivity computation
        //  As a Switch with delayed evaluation since adjointBoundaryConditions
        //  have not been allocated at the time of construction
        Switch hasBCdxdbMult_;


    // Protected Member Functions

        //- Compute, if necessary, and return the hasBCdxdbMult_ bool
        bool hasBCdxdbMult(const labelHashSet& sensitivityPatchIDs);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("incompressible");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            incompressibleAdjointSolver,
            dictionary,
            (
                fvMesh& mesh,
                const word& managerType,
                const dictionary& dict,
                const word& primalSolverName,
                const word& solverName
            ),
            (mesh, managerType, dict, primalSolverName, solverName)
        );


    // Constructors

        //- Construct from mesh and dictionary
        incompressibleAdjointSolver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    // Selectors

        //- Return a reference to the selected incompressible adjoint solver
        static autoPtr<incompressibleAdjointSolver> New
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    //- Destructor
    virtual ~incompressibleAdjointSolver() = default;


    // Member Functions

        // Access

            //- Access to the incompressible primal variables set
            const incompressibleVars& getPrimalVars() const;

            //- Access to the incompressible adjoint variables set
            virtual const incompressibleAdjointVars& getAdjointVars() const;

            //- Access to the incompressible adjoint variables set
            virtual incompressibleAdjointVars& getAdjointVars();

            //- Access to the ATC model
            const autoPtr<ATCModel>& getATCModel() const;

            //- Access to the ATC model
            autoPtr<ATCModel>& getATCModel();

            //- Should the adjoint to the eikonal equation be solved
            virtual bool includeDistance() const;

            //- Return the dimensions of the adjoint distance field
            virtual dimensionSet daDimensions() const;

            //- Return the dimensions of the adjoint grid displacement variable
            virtual dimensionSet maDimensions() const;

            //- Return the source the adjoint eikonal equation
            virtual tmp<volScalarField> adjointEikonalSource();

            //- Return the distance field, to be used in the solution of the
            //- adjoint eikonal PDE
            virtual tmp<volScalarField> yWall() const;


        // Evolution

            //- Update primal based quantities, e.g. the primal fields
            //- in adjoint turbulence models
            virtual void updatePrimalBasedQuantities();


        // Sensitivity related functions

            // Functions related to the computation of sensitivity derivatives
            // All functions get the field to accumulate their contibution on
            // as an argument. Should be implemented by the derived classes
            // if the physics there adds terms to the sensitivity derivatives

                // Shape optimisation

                    //- Compute the multiplier for grad(dxdb)
                    //  Used in shape sensitivity derivatives, computed with
                    //  the FI and E-SI approaches
                    virtual void accumulateGradDxDbMultiplier
                    (
                        volTensorField& gradDxDbMult,
                        const scalar dt
                    );

                    //- Compute the multiplier for div(dxdb)
                    //  Used in shape sensitivity derivatives, computed with
                    //  the FI and E-SI approaches
                    virtual void accumulateDivDxDbMultiplier
                    (
                        autoPtr<scalarField>& divDxDbMult,
                        const scalar dt
                    );

                    //- Accumulate the multipliers of geometric quantities
                    //- defined at the boundary, usually through an objective
                    //- or constraint function
                    virtual void accumulateGeometryVariationsMultipliers
                    (
                        autoPtr<boundaryVectorField>& dSfdbMult,
                        autoPtr<boundaryVectorField>& dnfdbMult,
                        autoPtr<boundaryVectorField>& dxdbDirectMult,
                        autoPtr<pointBoundaryVectorField>& pointDxDirectDbMult,
                        const labelHashSet& sensitivityPatchIDs,
                        const scalar dt
                    );

                    //- Contributions from boundary functions that inlcude
                    //- geometric aspects in them and change when the geometry
                    //- is displaced, e.g. rotationWallVelocity
                    virtual void accumulateBCSensitivityIntegrand
                    (
                        autoPtr<boundaryVectorField>& bcDxDbMult,
                        const labelHashSet& sensitivityPatchIDs,
                        const scalar dt
                    );

                    //- Contributions from fvOptions that inlcude
                    //- geometric aspects in them and change when the geometry
                    //- is displaced, e.g. MRF
                    virtual void accumulateOptionsDxDbMultiplier
                    (
                        vectorField& optionsDxDbMult,
                        const scalar dt
                    );


                // Topology optimisation

                    //- Compute the multiplier of beta
                    virtual void topOSensMultiplier
                    (
                        scalarField& betaMult,
                        const word& designVariablesName,
                        const scalar dt
                    );


        // IO

            //- In case of multi-point runs with turbulent flows,
            //- output dummy turbulence fields with the base names, to allow
            //- continuation
            virtual bool write(const bool valid = true) const
            {
                if (mesh_.time().writeTime())
                {
                    return primalVars_.write();
                }

                return false;
            }

            //- In case of multi-point runs with turbulent flows,
            //- output dummy turbulence fields with the base names, to allow
            //- continuation
            virtual bool writeNow() const
            {
                return primalVars_.write();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

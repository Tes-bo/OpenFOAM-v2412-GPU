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
    Foam::incompressiblePrimalSolver

Description
    Base class for primal incompressible solvers

\*---------------------------------------------------------------------------*/

#ifndef incompressiblePrimalSolver_H
#define incompressiblePrimalSolver_H

#include "primalSolver.H"
#include "incompressibleVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class objective;

/*---------------------------------------------------------------------------*\
                  Class incompressiblePrimalSolver Declaration
\*---------------------------------------------------------------------------*/

class incompressiblePrimalSolver
:
    public primalSolver
{
private:

    // Privare Member Functions

        //- No copy construct
        incompressiblePrimalSolver(const incompressiblePrimalSolver&) = delete;

        //- No copy assignment
        void operator=(const incompressiblePrimalSolver&) = delete;


protected:

    // Protected data

        //- Convergence criterion for reconstructing phi from U and p
        scalar phiReconstructionTol_;

        //- Max iterations  for reconstructing phi from U and p
        label phiReconstructionIters_;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("incompressible");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            incompressiblePrimalSolver,
            dictionary,
            (
                fvMesh& mesh,
                const word& managerType,
                const dictionary& dict,
                const word& solverName
            ),
            (mesh, managerType, dict, solverName)
        );


    // Constructors

        //- Construct from mesh and dictionary
        incompressiblePrimalSolver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& solverName
        );


    // Selectors

        //- Return a reference to the selected incompressible primal solver
        static autoPtr<incompressiblePrimalSolver> New
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& solverName
        );


    //- Destructor
    virtual ~incompressiblePrimalSolver() = default;


    // Member Functions

        //- Read dict if updated
        virtual bool readDict(const dictionary& dict);


        // Access

            //- Return the list of objectives assodicated with this solver
            UPtrList<objective> getObjectiveFunctions() const;

            //- Access to the incompressible variables set
            const incompressibleVars& getIncoVars() const;

            //- Access to the incompressible variables set
            incompressibleVars& getIncoVars();


        // Evolution

            //- Update boundary conditions
            virtual void correctBoundaryConditions();


        // IO

            //- In case of multi-point runs with turbulent flows,
            //- output dummy turbulence fields with the base names, to allow
            //- continuation
            virtual bool write(const bool valid = true) const
            {
                if (mesh_.time().writeTime())
                {
                    return getIncoVars().write();
                }

                return false;
            }

            //- In case of multi-point runs with turbulent flows,
            //- output dummy turbulence fields with the base names, to allow
            //- continuation
            virtual bool writeNow() const
            {
                return getIncoVars().write();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

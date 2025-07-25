/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::adjointMeshMovementSolver

Description
    Class solving the adjoint grid dispalcement PDEs.
    Assumes the primal grid displacement PDE is a Laplace one with uniform
    diffusivity.

    Reference:
    \verbatim
        For the derivation of the adjoint grid displacement PDEs, see
            Kavvadias, I., Papoutsis-Kiachagias, E., & Giannakoglou, K. (2015).
            On the proper treatment of grid sensitivities in continuous adjoint
            methods for shape optimization.
            Journal of Computational Physics, 301, 1–18.
            http://doi.org/10.1016/j.jcp.2015.08.012

    \endverbatim

SourceFiles
    adjointMeshMovementSolver.C

\*---------------------------------------------------------------------------*/

#ifndef adjointMeshMovementSolver_H
#define adjointMeshMovementSolver_H

#include "boundaryFieldsFwd.H"
#include "createZeroField.H"
#include "variablesSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class ShapeSensitivitiesBase;

/*---------------------------------------------------------------------------*\
                    Class adjointMeshMovementSolver Decleration
\*---------------------------------------------------------------------------*/

class adjointMeshMovementSolver
{
protected:

    // Protected Data Members

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Dictionary containing solution controls
        dictionary dict_;

        //- Part of sensitivity derivatives coming from the adjoint grid
        //- displacement PDE
        autoPtr<boundaryVectorField> meshMovementSensPtr_;

        // Underlaying adjoint sensitivities
        ShapeSensitivitiesBase& adjointSensitivity_;

        //- Adjoint grid displacement field
        volVectorField ma_;

        //- Source term of the adjoint grid displacement PDEs
        volVectorField source_;

        //- Solution controls
        label iters_;
        scalar tolerance_;


    // Protected Member Functions

        //- Read options each time a new solution is found
        void read();

        //- Set the source term of the PDE
        void setSource();


private:

    // Private Member Functions

        //- No copy construct
        adjointMeshMovementSolver(const adjointMeshMovementSolver&) = delete;

        //- No copy assignment
        void operator=( const adjointMeshMovementSolver) = delete;


public:

    //- Runtime type information
    TypeName("adjointMeshMovementSolver");


    // Constructors

        //- Construct from components
        adjointMeshMovementSolver
        (
            const fvMesh& mesh,
            const dictionary& dict,
            ShapeSensitivitiesBase& adjointSensitivity
        );


    // Destructor
    virtual ~adjointMeshMovementSolver() = default;


    // Member Functions

        //- Read dict if changed
        virtual bool readDict(const dictionary& dict);

        //- Calculate the adjoint distance field
        virtual void solve();

        //- Reset the source term
        void reset();

        //- Return the sensitivity term depending on ma
        boundaryVectorField& meshMovementSensitivities();

        //- Return the adjoint distance field
        inline const volVectorField& ma() const
        {
            return ma_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

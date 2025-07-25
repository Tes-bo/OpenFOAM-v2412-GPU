/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::adjointBoundaryCondition

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef adjointBoundaryCondition_H
#define adjointBoundaryCondition_H

#include "boundaryAdjointContribution.H"
#include "ATCModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class adjointBoundaryCondition Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class adjointBoundaryCondition
{
protected:

    // Protected data

        //- Reference to patch
        const fvPatch& patch_;

        //- objectiveManager name corresponding to field
        word managerName_;

        //- adjointSolver name corresponding to field
        word adjointSolverName_;

        //- simulationType corresponding to field.
        //  A placeholder for now
        word simulationType_;

        //- Engine to manage contributions of the objective functions
        //- to the adjoint boundary conditions
        autoPtr<boundaryAdjointContribution> boundaryContrPtr_;

        //- Whether to add the extra term from the UaGradU formulation.
        //  As a Switch with delayed evaluation since ATCModel
        //  has not been allocated at the time of construction
        Switch addATCUaGradUTerm_;


    // Protected Member Functions

        //- Get gradient of field on a specific boundary
        template<class Type2>
        tmp
        <
            Field<typename Foam::outerProduct<Foam::vector, Type2>::type>
        > computePatchGrad(word name);

        //- Whether to add the extra term from the UaGradU formulation
        bool addATCUaGradUTerm();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("adjointBoundaryCondition");


    // Constructors

        //- Construct from field and base name
        adjointBoundaryCondition
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF,
            const word& solverName
        );

        //- Construct as copy
        adjointBoundaryCondition(const adjointBoundaryCondition<Type>&);


    //- Destructor
    virtual ~adjointBoundaryCondition() = default;


    // Member Functions

        // Access

            //- Return objectiveManager name
            const word& objectiveManagerName() const;

            //- Return adjointSolverName
            const word& adjointSolverName() const;

            //- Return the simulationType
            const word& simulationType() const;

            //- Set the ptr to the correct boundaryAdjointContribution
            void setBoundaryContributionPtr();

            //- Get boundaryContribution
            boundaryAdjointContribution& getBoundaryAdjContribution();

            //- ATC type might be useful for a number of BCs. Return here
            const ATCModel& getATC() const;


        // Contribution to sensitivity derivatives

            //- Return contribution to sensitivity derivatives
            //  For adjoint boundary conditions corresponding to primal
            //  boundary conditions that include geometric components
            //  (e.g. rotatingWallVelocity)
            virtual tmp
            <
                Field<typename Foam::outerProduct<Foam::vector, Type>::type>
            > dxdbMult() const;

        //- Update the primal based quantities related to the adjoint boundary
        //- conditions
        virtual void updatePrimalBasedQuantities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "adjointBoundaryCondition.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

#define makeAdjointBoundaryCondition(adjointTypeBoundaryCondition)             \
defineNamedTemplateTypeNameAndDebug(adjointTypeBoundaryCondition, 0);

// ************************************************************************* //

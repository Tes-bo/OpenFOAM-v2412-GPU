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
    Foam::topODesignVariables

Description
    Design variables for porosity-based topology optimisation (topO) problems

    Reference:
    \verbatim
        The main idea behind porosity-based topO can be found in

            Papoutsis-Kiachagias, E. M., & Giannakoglou, K. C. (2014).
            Continuous Adjoint Methods for Turbulent Flows, Applied to Shape
            and Topology Optimization: Industrial Applications.
            Archives of Computational Methods in Engineering, 23(2), 255-299.
            http://doi.org/10.1007/s11831-014-9141-9

        with regularisation and projection approaches following

            Lazarov, B. S., & Sigmund O. (2010).
            Filters in topology optimization based on Helmholtz-type
            differential equations.
            International Journal for Numerical Methods in Engineering,
            86(6), 765-781.
            https://doi.org/10.1002/nme.3072
    \endverbatim

SourceFiles
    topODesignVariables.C

\*---------------------------------------------------------------------------*/

#ifndef topODesignVariables_H
#define topODesignVariables_H

#include "designVariables.H"
#include "topOVariablesBase.H"
#include "fieldRegularisation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class topODesignVariables Declaration
\*---------------------------------------------------------------------------*/

class topODesignVariables
:
    public topOVariablesBase,
    public designVariables
{
protected:

    // Protected data

        //- A subfield of the design variables corresponding to the porosity
        //- field.
        //  Usually the same as *this
        SubField<scalar> alpha_;

        //- Mechanism to regularise the field of design variables
        fieldRegularisation regularisation_;

        //- Write all fields related to imposition of the Brinkman
        //- penalisation (i.e. design variables, regularised design variables
        //- and the indicator field)
        bool writeAllFields_;

        //- Add the fvOptions necessary for topO automatically
        //  WIP
        bool addFvOptions_;


    // Protected Member Functions

        //- Update the design variables given their correction
        virtual void updateField
        (
            const scalarField& correction,
            const label fluidID = 0
        );

        //- Apply fixed values in certain zones
        void applyFixedValues();

        //- Set active design variables
        virtual void setActiveDesignVariables
        (
            const label fluidID = 0,
            const bool activeIO = false
        );

        //- Read field with (path of) the design variables and store input
        //- in the design variables list with optional offset
        void readField
        (
            const word& name,
            const label fluidID = 0,
            const bool setIOValues = false
        );

        //- Part of the constructor initialisation
        virtual void initialize();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topODesignVariables(const topODesignVariables&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const topODesignVariables&) = delete;


public:

    //- Runtime type information
    TypeName("topO");

    // Constructors

        //- Construct from dictionary
        topODesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct from dictionary and size
        topODesignVariables
        (
            fvMesh& mesh,
            const dictionary& dict,
            const label size
        );


    // Selectors

        //- Construct and return the selected design variables
        static autoPtr<topODesignVariables> New
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~topODesignVariables() = default;


    // Member Functions

        //- Get the indicator field
        virtual const volScalarField& beta() const;

        //- Return interpolant.
        //  Defaults to beta in mono-fluid topO optimization
        virtual const scalarField& interpolationField
        (
            const word& interpolationFieldName = "beta"
        ) const;

        //- Interpolate between the given field and solid values
        virtual void interpolate
        (
            volScalarField& field,
            const topOInterpolationFunction& interpolationFunc,
            const FieldField<Field, scalar>& fluidValues,
            const scalarField& solidValues,
            const label fieldi,
            const word& interpolationFieldName = "beta"
        ) const;

        //- Post-processing sensitivities due to interpolations based on
        //- the indicator fields
        virtual void interpolationSensitivities
        (
            scalarField& sens,
            const topOInterpolationFunction& interpolationFunc,
            const FieldField<Field, scalar>& fluidValues,
            const scalarField& solidValues,
            const label fieldi,
            const word& designVariablesName,
            const word& interpolationFieldName = "beta"
        ) const;

        //- Nullify given field in the solid domain
        virtual void nullifyInSolid
        (
            scalarField& field,
            const topOInterpolationFunction& interpolationFunc
        ) const;

        //- Nullify given field in the solid domain
        virtual void nullifyInSolidSensitivities
        (
            scalarField& sens,
            const topOInterpolationFunction& interpolationFunc,
            const word& designVariablesName
        ) const;

        //- Return the Brinkman penalisation term
        virtual tmp<scalarField> penalty
        (
            const word& interpolationFieldName,
            const topOInterpolationFunction& interpolationFunc
        ) const;

        //- Return the penalty term derivative
        virtual tmp<scalarField> penaltySensitivities
        (
            const word& interpolationFieldName,
            const topOInterpolationFunction& interpolationFunc
        ) const;

        //- Update design variables based on a given correction
        virtual void update(scalarField& correction);

        //- Compute eta if not set in the first step
        virtual scalar computeEta(scalarField& correction);

        //- Whether to use global sum when computing matrix-vector products
        //- in update methods
        virtual bool globalSum() const;

        //- Assemble sensitivity derivatives, by combining the part related
        //- to the primal and adjoint solution with the part related to the
        //- design variables.
        //  Takes into consideration the adjoint to the regularisation process
        //  and the derivative of the interpolationFunc
        virtual tmp<scalarField> assembleSensitivities
        (
            adjointSensitivity& sens
        );

        //- Set initial values of the design variables
        //  For design variables sets that need to be initialised after the
        //  construction of the primal fields.
        virtual void setInitialValues();

        //- Automatically add fvOptions depending on the design variables
        //- to the primal and adjoint solvers
        //  WIP
        virtual void addFvOptions
        (
            const PtrList<primalSolver>& primalSolver,
            const PtrList<adjointSolverManager>& adjointSolverManagers
        );

        //- Get access to the regularisation object
        inline fieldRegularisation& getRegularisation()
        {
            return regularisation_;
        }

        //- Write porosity field to file
        virtual void writeDesignVars();

        //- The writeData function required by the regIOobject write operation
        virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

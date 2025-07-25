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
    Foam::adjointOutletVelocityFluxFvPatchVectorField

Description
    An outlet boundary condition for patches in which the primal flow exhibits
    recirculation. Adds the contribution of the objective as an adjoint
    momentum flux directly to the PDEs, without the need to first compute an
    adjoint outlet velocity, circumventing thus the division with (almost) zero
    that manifests in case of primal flow recirculation.

SourceFiles
    adjointOutletVelocityFluxFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef adjointOutletVelocityFluxFvPatchVectorField_H
#define adjointOutletVelocityFluxFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "adjointBoundaryConditions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class adjointOutletVelocityFluxFvPatch Declaration
\*---------------------------------------------------------------------------*/

class adjointOutletVelocityFluxFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public adjointVectorBoundaryCondition
{

public:

    //- Runtime type information
    TypeName("adjointOutletVelocityFlux");


    // Constructors

        //- Construct from patch and internal field
        adjointOutletVelocityFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        adjointOutletVelocityFluxFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- adjointOutletVelocityFluxFvPatchVectorField
        //- onto a new patch
        adjointOutletVelocityFluxFvPatchVectorField
        (
            const adjointOutletVelocityFluxFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        adjointOutletVelocityFluxFvPatchVectorField
        (
            const adjointOutletVelocityFluxFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<vector>> clone() const
        {
            return fvPatchField<vector>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return fvPatchField<vector>::Clone(*this, iF);
        }


    // Member functions

        //- Return true: Allow adjoint solvers to obtain the outlet phia
        //  value through HbyA
        virtual bool assignable() const
        {
            return true;
        }

        //- add source term in the first cells off the wall due to adjoint WF
        virtual void manipulateMatrix(fvMatrix<vector>& matrix);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Return the matrix diagonal coefficients corresponding to the
        //- evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the
        //- evaluation of the value of this patchField with given weights
        virtual tmp<Field<vector>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&

        ) const;

        //- Return the matrix source coefficients corresponding to the
        //- evaluation of the gradient of this patchField
        virtual tmp<Field<vector>> gradientBoundaryCoeffs() const;

        //- Return the matrix diagonal coefficients corresponding to the
        //- evaluation of the gradient of this patchField
        virtual tmp<Field<vector>> gradientInternalCoeffs() const;

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::oversetFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    Boundary condition for use on overset patches. To be run in combination
    with special dynamicFvMesh type that inserts interpolation into the matrix.

SourceFiles
    oversetFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatchField_H
#define oversetFvPatchField_H

#include "oversetFvPatch.H"
#include "coupledFvPatchFields.H"
#include "oversetLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oversetFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class oversetFvPatchField
:
    public oversetLduInterfaceField,
    public coupledFvPatchField<Type>
{
protected:

    // Protected Data

        //- Local reference cast into the overset patch
        const oversetFvPatch& oversetPatch_;

        //- Master patch ID
        mutable label masterPatchID_;

        // Hole cell controls

            //- Flag to set hole cell values
            const bool setHoleCellValue_;

            //- Flag to correct fluxes
            const bool fluxCorrection_;

            //- Flag to interpolate hole cell values (from nearby non-hole cell)
            const bool interpolateHoleCellValue_;

            //- Hole cell value
            const Type holeCellValue_;

            //- Fringe upper coefficients
            mutable scalarField fringeUpperCoeffs_;

            //- Fringe lower coefficients
            mutable scalarField fringeLowerCoeffs_;

            //- Fringe faces
            mutable labelField fringeFaces_;

            //- Zone to sum flux for mass conservation
            label zoneId_;


public:

    //- Runtime type information
    TypeName(oversetFvPatch::typeName_());


    // Constructors

        //- Construct from patch and internal field
        oversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oversetFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given oversetFvPatchField onto a new patch
        oversetFvPatchField
        (
            const oversetFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oversetFvPatchField(const oversetFvPatchField<Type>&);

        //- Construct as copy setting internal field reference
        oversetFvPatchField
        (
            const oversetFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

        // Coupled and adjust flux

            //- Return neighbour field.  Dummy
            virtual tmp<Field<Type>> patchNeighbourField() const;

            //- Adjust psi for mass correction. Requires storeFringeCoefficients
            //  to have been called before
            void adjustPsi
            (
                solveScalarField& psi,
                const lduAddressing& lduAddr,
                solveScalarField& result
            ) const;

            //- Store fringe coefficients and faces
            void storeFringeCoefficients(const fvMatrix<Type>& matrix);

            //- Calculate patch flux (helper function). Requires
            //  storeFringeCoefficients to have been called before
            void fringeFlux
            (
                 const fvMatrix<Type>& m,
                 const surfaceScalarField& phi
            ) const;


        // Evaluation

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate(const Pstream::commsTypes commsType);

            //- Return the matrix diagonal coefficients corresponding to the
            //- evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type>>::New(this->size(), Foam::zero{});
            }

            //- Return the matrix source coefficients corresponding to the
            //- evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const
            {
                return tmp<Field<Type>>::New(this->size(), Foam::zero{});
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //- evaluation of the gradient of this patchField
            tmp<Field<Type>> gradientInternalCoeffs() const
            {
                return tmp<Field<Type>>::New(this->size(), Foam::zero{});
            }

            //- Return the matrix source coefficients corresponding to the
            //- evaluation of the gradient of this patchField
            tmp<Field<Type>> gradientBoundaryCoeffs() const
            {
                return tmp<Field<Type>>::New(this->size(), Foam::zero{});
            }

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<Type>& matrix);


        // Coupled interface functionality

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction,
                const Pstream::commsTypes commsType
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                Field<Type>&,
                const bool add,
                const lduAddressing& lduAddr,
                const label patchId,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                Field<Type>&,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const Field<Type>&,
                const scalarField&,
                const Pstream::commsTypes commsType
            ) const
            {
                NotImplemented;
            }

            virtual void initInterfaceMatrixUpdate
            (
                solveScalarField& result,
                const bool add,
                const lduAddressing&,
                const label interfacei,
                const solveScalarField& psiInternal,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "oversetFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

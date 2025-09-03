/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::fvMatrix

Description
    A special matrix type and solver, designed for finite volume
    solutions of scalar equations.
    Face addressing is used to make all matrix assembly
    and solution loops vectorise.

SourceFiles
    fvMatrix.C
    fvMatrixSolve.C
    fvScalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvMatrix_H
#define Foam_fvMatrix_H

#include "volFields.H"
#include "surfaceFields.H"
#include "lduMatrix.H"
#include "tmp.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "className.H"
#include "lduPrimitiveMeshAssembly.H"
#include "lduMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type> class fvMatrix;
template<class T> class UIndirectList;

template<class Type>
Ostream& operator<<(Ostream&, const fvMatrix<Type>&);


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> operator&
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> operator&
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> operator&
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> operator&
(
    const tmp<fvMatrix<Type>>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> operator&
(
    const tmp<fvMatrix<Type>>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> operator&
(
    const tmp<fvMatrix<Type>>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);


/*---------------------------------------------------------------------------*\
                           Class fvMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvMatrix
:
    public refCount,
    public lduMatrix
{
public:

    // Public Types

        //- The geometric field type for psi
        typedef
            GeometricField<Type, fvPatchField, volMesh>
            psiFieldType;

        //- Field type for face flux (for non-orthogonal correction)
        typedef
            GeometricField<Type, fvsPatchField, surfaceMesh>
            faceFluxFieldType;

        /// //- The internal field type for the psi field
        /// typedef DimensionedField<Type, volMesh> Internal;


private:

    // Private Data

        //- Const reference to field
        //  Converted into a non-const reference at the point of solution.
        const psiFieldType& psi_;

        //- Originating fvMatrices when assembling matrices. Empty if not used.
        PtrList<fvMatrix<Type>> subMatrices_;

        //- Is the fvMatrix using implicit formulation
        bool useImplicit_;

        //- Name of the lduAssembly
        word lduAssemblyName_;

        //- Number of fvMatrices added to this
        label nMatrix_;

        //- Dimension set
        dimensionSet dimensions_;

        //- Source term
        Field<Type> source_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //- for internal cells
        FieldField<Field, Type> internalCoeffs_;

        //- Boundary scalar field containing pseudo-matrix coeffs
        //- for boundary cells
        FieldField<Field, Type> boundaryCoeffs_;

        //- Face flux field for non-orthogonal correction
        mutable std::unique_ptr<faceFluxFieldType> faceFluxCorrectionPtr_;


protected:

    //- Declare friendship with the fvSolver class
    friend class fvSolver;


    // Protected Member Functions

        //- Add patch contribution to internal field
        template<class Type2>
        void addToInternalField
        (
            const labelUList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void addToInternalField
        (
            const labelUList& addr,
            const tmp<Field<Type2>>& tpf,
            Field<Type2>& intf
        ) const;

        //- Subtract patch contribution from internal field
        template<class Type2>
        void subtractFromInternalField
        (
            const labelUList& addr,
            const Field<Type2>& pf,
            Field<Type2>& intf
        ) const;

        template<class Type2>
        void subtractFromInternalField
        (
            const labelUList& addr,
            const tmp<Field<Type2>>& tpf,
            Field<Type2>& intf
        ) const;


        // Implicit helper functions

            //- Name the implicit assembly addressing
            //  \return true if assembly is used
            bool checkImplicit(const label fieldi = 0);


        // Matrix completion functionality

            void addBoundaryDiag
            (
                scalarField& diag,
                const direction cmpt
            ) const;

            void addCmptAvBoundaryDiag(scalarField& diag) const;

            void addBoundarySource
            (
                Field<Type>& source,
                const bool couples=true
            ) const;


        // Matrix manipulation functionality

            //- Set solution in given cells to the specified values
            template<template<class> class ListType>
            void setValuesFromList
            (
                const labelUList& cellLabels,
                const ListType<Type>& values
            );


public:

    //- Solver class returned by the solver function
    //- used for systems in which it is useful to cache the solver for reuse.
    //  E.g. if the solver is potentially expensive to construct (AMG) and can
    //  be used several times (PISO)
    class fvSolver
    {
        fvMatrix<Type>& fvMat_;

        autoPtr<lduMatrix::solver> solver_;

    public:

        // Constructors

            fvSolver(fvMatrix<Type>& fvMat, autoPtr<lduMatrix::solver>&& sol)
            :
                fvMat_(fvMat),
                solver_(std::move(sol))
            {}


        // Member Functions

            //- Solve returning the solution statistics.
            //  Solver controls read from dictionary
            SolverPerformance<Type> solve(const dictionary& solverControls);

            //- Solve returning the solution statistics.
            //  Solver controls read from fvSolution
            SolverPerformance<Type> solve();
    };


    // Runtime information
    ClassName("fvMatrix");


    // Constructors

        //- Construct given a field to solve for
        fvMatrix
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            const dimensionSet& ds
        );

        //- Copy construct
        fvMatrix(const fvMatrix<Type>&);

        //- Copy/move construct from tmp<fvMatrix<Type>>
        fvMatrix(const tmp<fvMatrix<Type>>&);

        //- Deprecated(2022-05) - construct with dimensionSet instead
        //  \deprecated(2022-05) - construct with dimensionSet instead
        FOAM_DEPRECATED_FOR(2022-05, "Construct with dimensionSet")
        fvMatrix
        (
            const GeometricField<Type, fvPatchField, volMesh>& psi,
            Istream& is
        )
        :
            fvMatrix<Type>(psi, dimensionSet(is))
        {}

        //- Construct and return a clone
        tmp<fvMatrix<Type>> clone() const
        {
            return tmp<fvMatrix<Type>>::New(*this);
        }


    //- Destructor
    virtual ~fvMatrix();


    // Member Functions

        // Access

            // Coupling

                label nMatrices() const
                {
                    return (nMatrix_ == 0 ? 1 : nMatrix_);
                }

                const fvMatrix<Type>& matrix(const label i) const
                {
                    return (nMatrix_ == 0 ? *this : subMatrices_[i]);
                }

                fvMatrix<Type>& matrix(const label i)
                {
                    return (nMatrix_ == 0 ? *this : subMatrices_[i]);
                }

                label globalPatchID
                (
                    const label fieldi,
                    const label patchi
                ) const
                {
                    if (!lduMeshPtr())
                    {
                        return patchi;
                    }
                    else
                    {
                        return lduMeshPtr()->patchMap()[fieldi][patchi];
                    }
                }

                //- Transfer lower, upper, diag and source to this fvMatrix
                void transferFvMatrixCoeffs();

                //- Create or update ldu assembly
                void createOrUpdateLduPrimitiveAssembly();

                //- Access to lduPrimitiveMeshAssembly
                lduPrimitiveMeshAssembly* lduMeshPtr();

                //- Const Access to lduPrimitiveMeshAssembly
                const lduPrimitiveMeshAssembly* lduMeshPtr() const;

                //- Manipulate matrix
                void manipulateMatrix(direction cmp);

                //- Manipulate boundary/internal coeffs for coupling
                void setBounAndInterCoeffs();

                //- Set interfaces
                void setInterfaces
                (
                    lduInterfaceFieldPtrsList&,
                    PtrDynList<lduInterfaceField>& newInterfaces
                );

                //- Add internal and boundary contribution to local patches
                void mapContributions
                (
                    label fieldi,
                    const FieldField<Field, Type>& fluxContrib,
                    FieldField<Field, Type>& contrib,
                    bool internal
                ) const;

                //- Return optional lduAdressing
                const lduPrimitiveMeshAssembly& lduMeshAssembly()
                {
                    return *lduMeshPtr();
                }

                //- Return psi
                const GeometricField<Type, fvPatchField, volMesh>& psi
                (
                    const label i = 0
                ) const
                {
                    return
                    (
                        (i == 0 && nMatrix_ == 0) ? psi_ : matrix(i).psi()
                    );
                }


                GeometricField<Type, fvPatchField, volMesh>& psi
                (
                    const label i = 0
                )
                {
                    return
                    (
                        (i == 0 && nMatrix_ == 0)
                      ? const_cast
                        <
                            GeometricField<Type, fvPatchField, volMesh>&
                        >(psi_)
                      : const_cast
                        <
                            GeometricField<Type, fvPatchField, volMesh>&
                        >
                        (
                            matrix(i).psi()
                        )
                    );
                }

                //- Clear multiple fvMatrices
                void clear()
                {
                    subMatrices_.clear();
                    nMatrix_ = 0;
                }


            const dimensionSet& dimensions() const noexcept
            {
                return dimensions_;
            }

            Field<Type>& source() noexcept
            {
                return source_;
            }

            const Field<Type>& source() const noexcept
            {
                return source_;
            }

            //- fvBoundary scalar field containing pseudo-matrix coeffs
            //- for internal cells
            const FieldField<Field, Type>& internalCoeffs() const noexcept
            {
                return internalCoeffs_;
            }

            //- fvBoundary scalar field containing pseudo-matrix coeffs
            //- for internal cells
            FieldField<Field, Type>& internalCoeffs() noexcept
            {
                return internalCoeffs_;
            }

            //- fvBoundary scalar field containing pseudo-matrix coeffs
            //- for boundary cells
            const FieldField<Field, Type>& boundaryCoeffs() const noexcept
            {
                return boundaryCoeffs_;
            }

            //- fvBoundary scalar field containing pseudo-matrix coeffs
            //- for boundary cells
            FieldField<Field, Type>& boundaryCoeffs() noexcept
            {
                return boundaryCoeffs_;
            }

            //- Declare return type of the faceFluxCorrectionPtr() function
            typedef std::unique_ptr<faceFluxFieldType> faceFluxFieldPtrType;

            //- Return pointer to face-flux non-orthogonal correction field
            faceFluxFieldPtrType& faceFluxCorrectionPtr()
            {
                return faceFluxCorrectionPtr_;
            }

            //- Set pointer to face-flux non-orthogonal correction field
            void faceFluxCorrectionPtr(faceFluxFieldType* flux)
            {
                faceFluxCorrectionPtr_.reset(flux);
            }

            //- True if face-flux non-orthogonal correction field exists
            bool hasFaceFluxCorrection() const noexcept
            {
                return bool(faceFluxCorrectionPtr_);
            }


        // Operations

            //- Set solution in given cells to the specified value
            //- and eliminate the corresponding equations from the matrix.
            void setValues
            (
                const labelUList& cellLabels,
                const Type& value
            );

            //- Set solution in given cells to the specified values
            //- and eliminate the corresponding equations from the matrix.
            void setValues
            (
                const labelUList& cellLabels,
                const UList<Type>& values
            );

            //- Set solution in given cells to the specified values
            //- and eliminate the corresponding equations from the matrix.
            void setValues
            (
                const labelUList& cellLabels,
                const UIndirectList<Type>& values
            );

            //- Set reference level for solution
            void setReference
            (
                const label celli,
                const Type& value,
                const bool forceReference = false
            );

            //- Set reference level for solution
            void setReferences
            (
                const labelUList& cellLabels,
                const Type& value,
                const bool forceReference = false
            );

            //- Set reference level for solution
            void setReferences
            (
                const labelUList& cellLabels,
                const UList<Type>& values,
                const bool forceReference = false
            );

            //- Set reference level for a component of the solution
            //- on a given patch face
            void setComponentReference
            (
                const label patchi,
                const label facei,
                const direction cmpt,
                const scalar value
            );

            //- Add fvMatrix
            void addFvMatrix(fvMatrix<Type>& matrix);

            //- Relax matrix (for steady-state solution).
            //  alpha = 1 : diagonally equal
            //  alpha < 1 : diagonally dominant
            //  alpha = 0 : do nothing
            //  Note: Requires positive diagonal.
            void relax(const scalar alpha);

            //- Relax matrix (for steady-state solution).
            //  alpha is read from controlDict
            void relax();

            //- Manipulate based on a boundary field
            void boundaryManipulate
            (
                typename GeometricField<Type, fvPatchField, volMesh>::
                    Boundary& values
            );

            //- Construct and return the solver
            //  Use the given solver controls
            autoPtr<fvSolver> solver(const dictionary&);

            //- Construct and return the solver
            //  Solver controls read from fvSolution
            autoPtr<fvSolver> solver();

            //- Solve segregated or coupled returning the solution statistics.
            //  Use the given solver controls
            SolverPerformance<Type> solveSegregatedOrCoupled(const dictionary&);

            //- Solve segregated or coupled returning the solution statistics.
            //  Solver controls read from fvSolution
            SolverPerformance<Type> solveSegregatedOrCoupled();

            //- Solve segregated returning the solution statistics.
            //  Use the given solver controls
            SolverPerformance<Type> solveSegregated(const dictionary&);

            //- Solve coupled returning the solution statistics.
            //  Use the given solver controls
            SolverPerformance<Type> solveCoupled(const dictionary&);

            //- Solve returning the solution statistics.
            //  Use the given solver controls
            SolverPerformance<Type> solve(const dictionary&);

            //- Solve returning the solution statistics.
            //  Uses \p name solver controls from fvSolution
            SolverPerformance<Type> solve(const word& name);

            //- Solve returning the solution statistics.
            //  Solver controls read from fvSolution
            SolverPerformance<Type> solve();

            //- Return the matrix residual
            tmp<Field<Type>> residual() const;

            //- Return the matrix scalar diagonal
            tmp<scalarField> D() const;

            //- Return the matrix Type diagonal
            tmp<Field<Type>> DD() const;

            //- Return the central coefficient
            tmp<volScalarField> A() const;

            //- Return the H operation source
            tmp<GeometricField<Type, fvPatchField, volMesh>> H() const;

            //- Return H(1)
            tmp<volScalarField> H1() const;

            //- Return the face-flux field from the matrix
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
                flux() const;

            //- Return the solver dictionary (from fvSolution) for \p name
            const dictionary& solverDict(const word& name) const;

            //- Return the solver dictionary for psi,
            //- taking into account finalIteration
            const dictionary& solverDict() const;


    // Member Operators

        void operator=(const fvMatrix<Type>&);
        void operator=(const tmp<fvMatrix<Type>>&);

        //- Inplace negate
        void negate();

        void operator+=(const fvMatrix<Type>&);
        void operator+=(const tmp<fvMatrix<Type>>&);

        void operator-=(const fvMatrix<Type>&);
        void operator-=(const tmp<fvMatrix<Type>>&);

        void operator+=(const DimensionedField<Type, volMesh>&);
        void operator+=(const tmp<DimensionedField<Type, volMesh>>&);
        void operator+=
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        );

        void operator-=(const DimensionedField<Type, volMesh>&);
        void operator-=(const tmp<DimensionedField<Type, volMesh>>&);
        void operator-=
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        );

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator+=(const Foam::zero) {}
        void operator-=(const Foam::zero) {}

        void operator*=(const volScalarField::Internal&);
        void operator*=(const tmp<volScalarField::Internal>&);
        void operator*=(const tmp<volScalarField>&);

        void operator*=(const dimensioned<scalar>&);


    // Friend Operators

        friend tmp<GeometricField<Type, fvPatchField, volMesh>>
        operator& <Type>
        (
            const fvMatrix<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        friend tmp<GeometricField<Type, fvPatchField, volMesh>>
        operator& <Type>
        (
            const fvMatrix<Type>&,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        );

        friend tmp<GeometricField<Type, fvPatchField, volMesh>>
        operator& <Type>
        (
            const tmp<fvMatrix<Type>>&,
            const DimensionedField<Type, volMesh>&
        );

        friend tmp<GeometricField<Type, fvPatchField, volMesh>>
        operator& <Type>
        (
            const tmp<fvMatrix<Type>>&,
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&
        );


    // Ostream Operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fvMatrix<Type>&
        );
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void checkMethod
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&,
    const char*
);

template<class Type>
void checkMethod
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&,
    const char*
);

template<class Type>
void checkMethod
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&,
    const char*
);


//- Solve returning the solution statistics given convergence tolerance
//  Use the given solver controls
template<class Type>
SolverPerformance<Type> solve
(
    fvMatrix<Type>&,
    const dictionary& solverControls
);

//- Solve returning the solution statistics given convergence tolerance,
//- deleting temporary matrix after solution.
//  Use the given solver controls
template<class Type>
SolverPerformance<Type> solve
(
    const tmp<fvMatrix<Type>>&,
    const dictionary& solverControls
);


//- Solve returning the solution statistics given convergence tolerance
//  Uses \p name solver controls from fvSolution
template<class Type>
SolverPerformance<Type> solve(fvMatrix<Type>&, const word& name);

//- Solve returning the solution statistics given convergence tolerance,
//- deleting temporary matrix after solution.
//  Uses \p name solver controls from fvSolution
template<class Type>
SolverPerformance<Type> solve(const tmp<fvMatrix<Type>>&, const word& name);


//- Solve returning the solution statistics given convergence tolerance
//  Uses solver controls from fvSolution
template<class Type>
SolverPerformance<Type> solve(fvMatrix<Type>&);

//- Solve returning the solution statistics given convergence tolerance,
//- deleting temporary matrix after solution.
//  Uses solver controls from fvSolution
template<class Type>
SolverPerformance<Type> solve(const tmp<fvMatrix<Type>>&);


//- Return the correction form of the given matrix
//- by subtracting the matrix multiplied by the current field
template<class Type>
tmp<fvMatrix<Type>> correction(const fvMatrix<Type>&);


//- Return the correction form of the given temporary matrix
//- by subtracting the matrix multiplied by the current field
template<class Type>
tmp<fvMatrix<Type>> correction(const tmp<fvMatrix<Type>>&);


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const tmp<fvMatrix<Type>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const fvMatrix<Type>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const tmp<fvMatrix<Type>>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const tmp<fvMatrix<Type>>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const tmp<fvMatrix<Type>>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const tmp<fvMatrix<Type>>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const tmp<fvMatrix<Type>>&,
    const dimensioned<Type>&
);


template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const fvMatrix<Type>&,
    const Foam::zero
);

template<class Type>
tmp<fvMatrix<Type>> operator==
(
    const tmp<fvMatrix<Type>>&,
    const Foam::zero
);


//- Unary negation
template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const fvMatrix<Type>&
);

//- Unary negation
template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<fvMatrix<Type>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const fvMatrix<Type>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<fvMatrix<Type>>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<fvMatrix<Type>>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<fvMatrix<Type>>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<fvMatrix<Type>>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const DimensionedField<Type, volMesh>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const DimensionedField<Type, volMesh>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const tmp<fvMatrix<Type>>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const dimensioned<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator+
(
    const dimensioned<Type>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const fvMatrix<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const fvMatrix<Type>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const fvMatrix<Type>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const fvMatrix<Type>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const fvMatrix<Type>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&,
    const DimensionedField<Type, volMesh>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&,
    const tmp<DimensionedField<Type, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const DimensionedField<Type, volMesh>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const DimensionedField<Type, volMesh>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<DimensionedField<Type, volMesh>>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const fvMatrix<Type>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const tmp<fvMatrix<Type>>&,
    const dimensioned<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const dimensioned<Type>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator-
(
    const dimensioned<Type>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const volScalarField::Internal&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const tmp<volScalarField::Internal>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const tmp<volScalarField>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const volScalarField::Internal&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const tmp<volScalarField::Internal>&,
    const tmp<fvMatrix<Type>>&
);

template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const tmp<volScalarField>&,
    const tmp<fvMatrix<Type>>&
);


template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const dimensioned<scalar>&,
    const fvMatrix<Type>&
);

template<class Type>
tmp<fvMatrix<Type>> operator*
(
    const dimensioned<scalar>&,
    const tmp<fvMatrix<Type>>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMatrix.C"
#endif

// Specialisation for scalars
#include "fvScalarMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

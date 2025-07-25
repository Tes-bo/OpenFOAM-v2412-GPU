/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
    Foam::LduMatrix

Description
    LduMatrix is a general matrix class in which the coefficients are
    stored as three arrays, one for the upper triangle, one for the
    lower triangle and a third for the diagonal.

    Addressing arrays must be supplied for the upper and lower triangles.

Note
    It might be better if this class were organised as a hierarchy starting
    from an empty matrix, then deriving diagonal, symmetric and asymmetric
    matrices.

SourceFiles
    LduMatrixATmul.C
    LduMatrix.C
    LduMatrixOperations.C
    LduMatrixSolver.C
    LduMatrixPreconditioner.C
    LduMatrixTests.C
    LduMatrixUpdateMatrixInterfaces.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LduMatrix_H
#define Foam_LduMatrix_H

#include "lduMesh.H"
#include "lduMatrix.H"
#include "Field.H"
#include "FieldField.H"
#include "LduInterfaceFieldPtrsList.H"
#include "SolverPerformance.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Type, class DType, class LUType> class LduMatrix;

template<class Type, class DType, class LUType>
Ostream& operator<<
(
    Ostream&,
    const LduMatrix<Type, DType, LUType>&
);


/*---------------------------------------------------------------------------*\
                           Class LduMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class DType, class LUType>
class LduMatrix
{
    // Private Data

        //- LDU mesh reference
        const lduMesh& lduMesh_;

        //- Diagonal coefficients
        std::unique_ptr<Field<DType>> diagPtr_;

        //- Off-diagonal coefficients
        std::unique_ptr<Field<LUType>> upperPtr_;

        //- Off-diagonal coefficients
        std::unique_ptr<Field<LUType>> lowerPtr_;

        //- Source
        std::unique_ptr<Field<Type>> sourcePtr_;

        //- Field interfaces (processor patches etc.)
        LduInterfaceFieldPtrsList<Type> interfaces_;

        //- Off-diagonal coefficients for interfaces
        FieldField<Field, LUType> interfacesUpper_;

        //- Off-diagonal coefficients for interfaces
        FieldField<Field, LUType> interfacesLower_;


public:

    friend class SolverPerformance<Type>;

    // -----------------------------------------------------------------------
    //- Abstract base-class for LduMatrix solvers
    class solver
    {
    protected:

        // Protected Data

            word fieldName_;
            const LduMatrix<Type, DType, LUType>& matrix_;

            //- Dictionary of solution controls
            dictionary controlDict_;

            //- Verbosity level for solver output statements
            int log_;

            //- Minimum number of iterations in the solver
            label minIter_;

            //- Maximum number of iterations in the solver
            label maxIter_;

            //- The matrix normalisation type
            lduMatrix::normTypes normType_;

            //- Final convergence tolerance
            Type tolerance_;

            //- Convergence tolerance relative to the initial
            Type relTol_;


        // Protected Member Functions

            //- Read the control parameters from controlDict_
            virtual void readControls();


        // Housekeeping

            //- Deprecated(2021-09) Read control parameter from dictionary
            //  \deprecated(2021-09) - use dictionary methods directly
            template<class T>
            void readControl(const dictionary& dict, T& val, const word& key)
            {
                dict.readIfPresent(key, val);
            }


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                symMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix,
                    const dictionary& solverDict
                ),
                (
                    fieldName,
                    matrix,
                    solverDict
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                solver,
                asymMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix,
                    const dictionary& solverDict
                ),
                (
                    fieldName,
                    matrix,
                    solverDict
                )
            );


        // Constructors

            //- Construct for given field name, matrix and controls
            solver
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix,
                const dictionary& solverDict
            );


        // Selectors

            //- Return a new solver
            static autoPtr<solver> New
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix,
                const dictionary& solverDict
            );


        //- Destructor
        virtual ~solver() = default;


        // Member Functions

            const word& fieldName() const noexcept
            {
                return fieldName_;
            }

            const LduMatrix<Type, DType, LUType>& matrix() const noexcept
            {
                return matrix_;
            }


            //- Read and reset the solver parameters from the given dictionary
            virtual void read(const dictionary&);

            virtual SolverPerformance<Type> solve
            (
                Field<Type>& psi
            ) const = 0;

            //- Return the matrix norm using the specified norm method
            Type normFactor
            (
                const Field<Type>& psi,
                const Field<Type>& Apsi,
                Field<Type>& tmpField,
                const lduMatrix::normTypes normType
            ) const;

            //- Return the matrix norm used to normalise the residual for the
            //- stopping criterion
            Type normFactor
            (
                const Field<Type>& psi,
                const Field<Type>& Apsi,
                Field<Type>& tmpField
            ) const
            {
                return this->normFactor(psi, Apsi, tmpField, normType_);
            }
    };


    // -----------------------------------------------------------------------
    //- Abstract base-class for LduMatrix smoothers
    class smoother
    {
    protected:

        // Protected Data

            word fieldName_;
            const LduMatrix<Type, DType, LUType>& matrix_;


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                symMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix
                ),
                (
                    fieldName,
                    matrix
                )
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                smoother,
                asymMatrix,
                (
                    const word& fieldName,
                    const LduMatrix<Type, DType, LUType>& matrix
                ),
                (
                    fieldName,
                    matrix
                )
            );


        // Constructors

            //- Construct for given field name and matrix
            smoother
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix
            );


        // Selectors

            //- Return a new smoother
            static autoPtr<smoother> New
            (
                const word& fieldName,
                const LduMatrix<Type, DType, LUType>& matrix,
                const dictionary& smootherDict
            );


        //- Destructor
        virtual ~smoother() = default;


        // Member Functions

            const word& fieldName() const noexcept
            {
                return fieldName_;
            }

            const LduMatrix<Type, DType, LUType>& matrix() const noexcept
            {
                return matrix_;
            }


            //- Smooth the solution for a given number of sweeps
            virtual void smooth
            (
                Field<Type>& psi,
                const label nSweeps
            ) const = 0;
    };


    // -----------------------------------------------------------------------
    //- Abstract base-class for LduMatrix preconditioners
    class preconditioner
    {
    protected:

        // Protected Data

            //- Reference to the base-solver this preconditioner is used with
            const solver& solver_;


    public:

        //- Runtime type information
        virtual const word& type() const = 0;


        // Declare run-time constructor selection tables

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                symMatrix,
                (
                    const solver& sol,
                    const dictionary& preconditionerDict
                ),
                (sol, preconditionerDict)
            );

            declareRunTimeSelectionTable
            (
                autoPtr,
                preconditioner,
                asymMatrix,
                (
                    const solver& sol,
                    const dictionary& preconditionerDict
                ),
                (sol, preconditionerDict)
            );


        // Constructors

            //- Construct for given solver
            preconditioner(const solver& sol)
            :
                solver_(sol)
            {}


        // Selectors

            //- Return a new preconditioner
            static autoPtr<preconditioner> New
            (
                const solver& sol,
                const dictionary& preconditionerDict
            );


        //- Destructor
        virtual ~preconditioner() = default;


        // Member functions

            //- Read and reset the preconditioner parameters
            //- from the given dictionary
            virtual void read(const dictionary&)
            {}

            //- Return wA the preconditioned form of residual rA
            virtual void precondition
            (
                Field<Type>& wA,
                const Field<Type>& rA
            ) const = 0;

            //- Return wT the transpose-matrix preconditioned form of
            //- residual rT.
            //  This is only required for preconditioning asymmetric matrices.
            virtual void preconditionT
            (
                Field<Type>& wT,
                const Field<Type>& rT
            ) const
            {
                NotImplemented;
            }
    };


    // -----------------------------------------------------------------------

    // Static Data

        // Declare name of the class and its debug switch
        ClassName("LduMatrix");


    // Constructors

        //- Construct given an LDU addressed mesh.
        //  The coefficients are initially empty for subsequent setting.
        //  Not yet 'explicit' (legacy code may rely on implicit construct)
        LduMatrix(const lduMesh& mesh);

        //- Copy construct
        LduMatrix(const LduMatrix<Type, DType, LUType>&);

        //- Move construct
        LduMatrix(LduMatrix<Type, DType, LUType>&&);

        //- Construct as copy or re-use as specified.
        LduMatrix(LduMatrix<Type, DType, LUType>&, bool reuse);

        //- Construct given an LDU addressed mesh and an Istream
        //- from which the coefficients are read
        LduMatrix(const lduMesh& mesh, Istream& is);


    //- Destructor
    ~LduMatrix() = default;


    // Member Functions

    // Addressing

        //- Return the LDU mesh from which the addressing is obtained
        const lduMesh& mesh() const noexcept
        {
            return lduMesh_;
        }

        //- Return the LDU addressing
        const lduAddressing& lduAddr() const
        {
            return lduMesh_.lduAddr();
        }

        //- Return the patch evaluation schedule
        const lduSchedule& patchSchedule() const
        {
            return lduMesh_.lduAddr().patchSchedule();
        }

        //- Const access to the interfaces
        const LduInterfaceFieldPtrsList<Type>& interfaces() const noexcept
        {
            return interfaces_;
        }

        //- Non-const access to the interfaces
        LduInterfaceFieldPtrsList<Type>& interfaces() noexcept
        {
            return interfaces_;
        }


    // Coefficients

        const Field<DType>& diag() const;
        const Field<LUType>& upper() const;
        const Field<LUType>& lower() const;
        const Field<Type>& source() const;

        Field<DType>& diag();
        Field<LUType>& upper();
        Field<LUType>& lower();
        Field<Type>& source();


    // Interfaces

        const FieldField<Field, LUType>& interfacesUpper() const noexcept
        {
            return interfacesUpper_;
        }

        const FieldField<Field, LUType>& interfacesLower() const noexcept
        {
            return interfacesLower_;
        }

        FieldField<Field, LUType>& interfacesUpper() noexcept
        {
            return interfacesUpper_;
        }

        FieldField<Field, LUType>& interfacesLower() noexcept
        {
            return interfacesLower_;
        }


    // Characteristics

        //- The matrix type (empty, diagonal, symmetric, ...)
        word matrixTypeName() const;

        bool hasDiag() const noexcept { return bool(diagPtr_); }
        bool hasUpper() const noexcept { return bool(upperPtr_); }
        bool hasLower() const noexcept { return bool(lowerPtr_); }
        bool hasSource() const noexcept { return bool(sourcePtr_); }

        //- Matrix has diagonal only
        bool diagonal() const noexcept
        {
            return (diagPtr_ && !lowerPtr_ && !upperPtr_);
        }

        //- Matrix is symmetric
        bool symmetric() const noexcept
        {
            return (diagPtr_ && !lowerPtr_ && upperPtr_);
        }

        //- Matrix is asymmetric (ie, full)
        bool asymmetric() const noexcept
        {
            return (diagPtr_ && lowerPtr_ && upperPtr_);
        }


        // Operations

            void sumDiag();
            void negSumDiag();

            void sumMagOffDiag(Field<LUType>& sumOff) const;

            //- Matrix multiplication
            void Amul(Field<Type>&, const tmp<Field<Type>>&) const;

            //- Matrix transpose multiplication
            void Tmul(Field<Type>&, const tmp<Field<Type>>&) const;


            //- Sum the coefficients on each row of the matrix
            void sumA(Field<Type>&) const;


            void residual(Field<Type>& rA, const Field<Type>& psi) const;

            tmp<Field<Type>> residual(const Field<Type>& psi) const;


            //- Initialise the update of interfaced interfaces
            //  for matrix operations
            void initMatrixInterfaces
            (
                const bool add,
                const FieldField<Field, LUType>& interfaceCoeffs,
                const Field<Type>& psiif,
                Field<Type>& result
            ) const;

            //- Update interfaced interfaces for matrix operations
            void updateMatrixInterfaces
            (
                const bool add,
                const FieldField<Field, LUType>& interfaceCoeffs,
                const Field<Type>& psiif,
                Field<Type>& result,
                const label startRequest // starting request (for non-blocking)
            ) const;


            tmp<Field<Type>> H(const Field<Type>&) const;
            tmp<Field<Type>> H(const tmp<Field<Type>>&) const;

            tmp<Field<Type>> faceH(const Field<Type>&) const;
            tmp<Field<Type>> faceH(const tmp<Field<Type>>&) const;


    // Member Operators

        //- Copy assignment
        void operator=(const LduMatrix<Type, DType, LUType>&);

        //- Move assignment
        void operator=(LduMatrix<Type, DType, LUType>&&);

        void negate();

        void operator+=(const LduMatrix<Type, DType, LUType>&);
        void operator-=(const LduMatrix<Type, DType, LUType>&);

        void operator*=(const scalarField&);
        void operator*=(scalar);


    // Ostream Operator

        friend Ostream& operator<< <Type, DType, LUType>
        (
            Ostream&,
            const LduMatrix<Type, DType, LUType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeLduMatrix(Type, DType, LUType)                                     \
                                                                               \
typedef Foam::LduMatrix<Type, DType, LUType>                                   \
    ldu##Type##DType##LUType##Matrix;                                          \
                                                                               \
defineNamedTemplateTypeNameAndDebug(ldu##Type##DType##LUType##Matrix, 0);      \
                                                                               \
                                                                               \
typedef LduMatrix<Type, DType, LUType>::smoother                               \
    ldu##Type##DType##LUType##Smoother;                                        \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Smoother,                                        \
    symMatrix                                                                  \
);                                                                             \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Smoother,                                        \
    asymMatrix                                                                 \
);                                                                             \
                                                                               \
                                                                               \
typedef LduMatrix<Type, DType, LUType>::preconditioner                         \
    ldu##Type##DType##LUType##Preconditioner;                                  \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Preconditioner,                                  \
    symMatrix                                                                  \
);                                                                             \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Preconditioner,                                  \
    asymMatrix                                                                 \
);                                                                             \
                                                                               \
                                                                               \
typedef LduMatrix<Type, DType, LUType>::solver                                 \
    ldu##Type##DType##LUType##Solver;                                          \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Solver,                                          \
    symMatrix                                                                  \
);                                                                             \
                                                                               \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ldu##Type##DType##LUType##Solver,                                          \
    asymMatrix                                                                 \
);


#define makeLduPreconditioner(Precon, Type, DType, LUType)                     \
                                                                               \
typedef Precon<Type, DType, LUType>                                            \
    Precon##Type##DType##LUType##Preconditioner;                               \
defineNamedTemplateTypeNameAndDebug                                            \
(                                                                              \
    Precon##Type##DType##LUType##Preconditioner,                               \
    0                                                                          \
);

#define makeLduSymPreconditioner(Precon, Type, DType, LUType)                  \
                                                                               \
LduMatrix<Type, DType, LUType>::preconditioner::                               \
addsymMatrixConstructorToTable<Precon##Type##DType##LUType##Preconditioner>    \
add##Precon##Type##DType##LUType##PreconditionerSymMatrixConstructorToTable_;

#define makeLduAsymPreconditioner(Precon, Type, DType, LUType)                 \
                                                                               \
LduMatrix<Type, DType, LUType>::preconditioner::                               \
addasymMatrixConstructorToTable<Precon##Type##DType##LUType##Preconditioner>   \
add##Precon##Type##DType##LUType##PreconditionerAsymMatrixConstructorToTable_;


#define makeLduSmoother(Smoother, Type, DType, LUType)                         \
                                                                               \
typedef Smoother<Type, DType, LUType>                                          \
    Smoother##Type##DType##LUType##Smoother;                                   \
                                                                               \
defineNamedTemplateTypeNameAndDebug                                            \
(                                                                              \
    Smoother##Type##DType##LUType##Smoother,                                   \
    0                                                                          \
);

#define makeLduSymSmoother(Smoother, Type, DType, LUType)                      \
                                                                               \
LduMatrix<Type, DType, LUType>::smoother::                                     \
    addsymMatrixConstructorToTable<Smoother##Type##DType##LUType##Smoother>    \
    add##Smoother##Type##DType##LUType##SymMatrixConstructorToTable_;

#define makeLduAsymSmoother(Smoother, Type, DType, LUType)                     \
                                                                               \
LduMatrix<Type, DType, LUType>::smoother::                                     \
    addasymMatrixConstructorToTable<Smoother##Type##DType##LUType##Smoother>   \
    add##Smoother##Type##DType##LUType##AsymMatrixConstructorToTable_;


#define makeLduSolver(Solver, Type, DType, LUType)                             \
                                                                               \
typedef Solver<Type, DType, LUType>                                            \
    Solver##Type##DType##LUType##Solver;                                       \
                                                                               \
defineNamedTemplateTypeNameAndDebug                                            \
(                                                                              \
    Solver##Type##DType##LUType##Solver,                                       \
    0                                                                          \
);

#define makeLduSymSolver(Solver, Type, DType, LUType)                          \
                                                                               \
LduMatrix<Type, DType, LUType>::solver::                                       \
    addsymMatrixConstructorToTable<Solver##Type##DType##LUType##Solver>        \
    add##Solver##Type##DType##LUType##SymMatrixConstructorToTable_;

#define makeLduAsymSolver(Solver, Type, DType, LUType)                         \
                                                                               \
LduMatrix<Type, DType, LUType>::solver::                                       \
    addasymMatrixConstructorToTable<Solver##Type##DType##LUType##Solver>       \
    add##Solver##Type##DType##LUType##AsymMatrixConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LduMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

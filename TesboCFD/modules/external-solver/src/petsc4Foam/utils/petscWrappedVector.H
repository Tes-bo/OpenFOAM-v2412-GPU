/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::PetscWrappedVector

Description
    Wrap an OpenFOAM list as a PETSc vector

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscWrappedVector_H
#define Foam_module_petscWrappedVector_H

// OpenFOAM
#include "List.H"
#include "scalar.H"
#include "vector.H"

// Interfacing
#include "petscErrorHandling.H"

// PETSc
#include "petscvec.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class PetscWrappedVector Declaration
\*---------------------------------------------------------------------------*/

class PetscWrappedVector
{
    // Private Data

        //- The PETSc vector
        Vec vec_;

        //- Call VecResetArray on destruction?
        bool reset_;


    // Private Member Functions

        //- Wrap an OpenFOAM list of values as a PETSc vector
        template<class T>
        void createImpl
        (
            const UList<T>& list,
            Vec vec
        )
        {
            PetscInt oldSize = 0;
            AssertPETSc(VecGetLocalSize(vec, &oldSize));
            if (list.size() != oldSize)
            {
                FatalErrorInFunction
                    << "List size " << list.size()
                    << " != vector size " << oldSize
                    << exit(FatalError);
            }
            AssertPETSc(PetscObjectReference((PetscObject)vec));
            vec_ = vec;
            AssertPETSc(VecPlaceArray
            (
                vec_,
                reinterpret_cast<solveScalar*>(const_cast<T*>(list.cdata()))
            ));
            reset_ = true;
        }

        //- Wrap an OpenFOAM list of values as a PETSc vector
        template<class T>
        void createImpl
        (
            const UList<T>& list,
            Mat mat
        )
        {
            AssertPETSc(MatCreateVecs(mat, &vec_, NULL));
            AssertPETSc(VecPlaceArray
            (
                vec_,
                reinterpret_cast<solveScalar*>(const_cast<T*>(list.cdata()))
            ));
            reset_ = true;
        }

        //- Wrap an OpenFOAM list of values as a PETSc vector
        template<class T>
        void createImpl
        (
            const UList<T>& list,
            MPI_Comm comm
        )
        {
            AssertPETSc(VecCreateMPIWithArray
            (
                comm,
                pTraits<T>::nComponents,
                list.size(),
                PETSC_DECIDE,
                reinterpret_cast<solveScalar*>(const_cast<T*>(list.cdata())),
                &vec_
            ));
            reset_ = false;
        }


    // Constructors

        //- Default construct
        PetscWrappedVector()
        :
            vec_(),
            reset_(false)
        {}


    // Generated Methods

        //- No copy construct
        PetscWrappedVector(const PetscWrappedVector&) = delete;

        //- No copy assignment
        void operator=(const PetscWrappedVector&) = delete;


public:

    // Constructors

        //- Wrap an OpenFOAM list of values as a PETSc vector
        //  Reuses an already created PETSc vector
        PetscWrappedVector
        (
            const UList<solveScalar>& list,
            Vec vec
        )
        :
            PetscWrappedVector()
        {
            createImpl(list, vec);
        }


        //- Wrap an OpenFOAM list of values as a PETSc vector
        //  Creates a new PETSc vector
        PetscWrappedVector
        (
            const UList<solveScalar>& list,
            Mat mat
        )
        :
            PetscWrappedVector()
        {
            createImpl(list, mat);
        }

        //- Wrap an OpenFOAM list of values as a PETSc vector
        explicit PetscWrappedVector
        (
            const UList<solveScalar>& list,
            MPI_Comm comm = PETSC_COMM_WORLD
        )
        :
            PetscWrappedVector()
        {
            createImpl(list, comm);
        }


    //- Destructor
    ~PetscWrappedVector()
    {
        if (reset_)
        {
            AssertPETSc(VecResetArray(vec_));
        }
        AssertPETSc(VecDestroy(&vec_));
    }


    // Member Operators

        //- Allow cast to a const PETSc vector
        operator const Vec&() const noexcept
        {
            return vec_;
        }

        //- Allow cast to a non-const PETSc vector
        operator Vec&() noexcept
        {
            return vec_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

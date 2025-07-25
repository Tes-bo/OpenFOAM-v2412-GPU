/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::SphericalTensor

Description
    A templated (3 x 3) diagonal tensor of objects of \<T\>, effectively
    containing 1 element, derived from VectorSpace.

See also
    Test-SphericalTensor.C

SourceFiles
    SphericalTensorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_SphericalTensor_H
#define Foam_SphericalTensor_H

#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SphericalTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SphericalTensor
:
    public VectorSpace<SphericalTensor<Cmpt>, Cmpt, 1>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef SphericalTensor<label> labelType;


    // Member Constants

        //- Rank of SphericalTensor is 2
        static constexpr direction rank = 2;


    // Static Data Members

        static const SphericalTensor I;
        static const SphericalTensor oneThirdI;
        static const SphericalTensor twoThirdsI;


    //- Component labeling enumeration
    enum components { II };


    // Generated Methods

        //- Default construct
        SphericalTensor() = default;

        //- Copy construct
        SphericalTensor(const SphericalTensor&) = default;

        //- Copy assignment
        SphericalTensor& operator=(const SphericalTensor&) = default;


    // Constructors

        //- Construct initialized to zero
        inline SphericalTensor(const Foam::zero);

        //- Construct given VectorSpace
        template<class Cmpt2>
        inline SphericalTensor
        (
            const VectorSpace<SphericalTensor<Cmpt2>, Cmpt2, 1>&
        );

        //- Construct given the component
        inline explicit SphericalTensor(const Cmpt& tii);

        //- Construct from Istream
        inline explicit SphericalTensor(Istream& is);


    // Member Functions

    // Component Access

        const Cmpt& ii() const noexcept { return this->v_[II]; }

        Cmpt& ii() noexcept { return this->v_[II]; }


    // Diagonal access and manipulation

        //- The L2-norm squared of the diagonal
        inline scalar diagSqr() const;


    // Tensor Operations

        //- Return non-Hermitian transpose (no-op)
        const SphericalTensor<Cmpt>& T() const noexcept { return *this; }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<SphericalTensor<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<SphericalTensor<Cmpt>>
:
    is_contiguous_label<Cmpt>
{};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<SphericalTensor<Cmpt>>
:
    is_contiguous_scalar<Cmpt>
{};


template<class Cmpt>
class typeOfSolve<SphericalTensor<Cmpt>>
{
public:

    typedef SphericalTensor<solveScalar> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SphericalTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

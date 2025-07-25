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
    Foam::SphericalTensor2D

Description
    A templated (2 x 2) diagonal tensor of objects of \<T\>, effectively
    containing 1 element, derived from VectorSpace.

See also
    Test-SphericalTensor2D.C

SourceFiles
    SphericalTensor2DI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_SphericalTensor2D_H
#define Foam_SphericalTensor2D_H

#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SphericalTensor2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SphericalTensor2D
:
    public VectorSpace<SphericalTensor2D<Cmpt>, Cmpt, 1>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef SphericalTensor2D<label> labelType;


    // Member Constants

        //- Rank of SphericalTensor2D is 2
        static constexpr direction rank = 2;


    // Static Data Members

        static const SphericalTensor2D I;
        static const SphericalTensor2D oneThirdI;
        static const SphericalTensor2D twoThirdsI;


    //- Component labeling enumeration
    enum components { II };


    // Generated Methods

        //- Default construct
        SphericalTensor2D() = default;

        //- Copy construct
        SphericalTensor2D(const SphericalTensor2D&) = default;

        //- Copy assignment
        SphericalTensor2D& operator=(const SphericalTensor2D&) = default;


    // Constructors

        //- Construct initialized to zero
        inline SphericalTensor2D(const Foam::zero);

        //- Construct given VectorSpace
        inline SphericalTensor2D
        (
            const VectorSpace<SphericalTensor2D<Cmpt>, Cmpt, 1>&
        );

        //- Construct given the component
        inline explicit SphericalTensor2D(const Cmpt& tii);

        //- Construct from Istream
        inline explicit SphericalTensor2D(Istream& is);


    // Member Functions

    // Component Access

        const Cmpt& ii() const noexcept { return this->v_[II]; }

        Cmpt& ii() noexcept { return this->v_[II]; }


    // Diagonal access and manipulation

        //- The L2-norm squared of the diagonal
        inline scalar diagSqr() const;


    // Tensor Operations

        //- Return non-Hermitian transpose (no-op)
        const SphericalTensor2D<Cmpt>& T() const noexcept { return *this; }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<SphericalTensor2D<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<SphericalTensor2D<Cmpt>>
:
    is_contiguous_label<Cmpt>
{};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<SphericalTensor2D<Cmpt>>
:
    is_contiguous_scalar<Cmpt>
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SphericalTensor2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

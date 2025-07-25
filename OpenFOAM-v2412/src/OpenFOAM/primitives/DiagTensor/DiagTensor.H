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
    Foam::DiagTensor

Description
    A templated (3 x 3) diagonal tensor of objects of \<T\>, effectively
    containing 3 elements, derived from VectorSpace.

See also
    Test-DiagTensor.C

SourceFiles
    DiagTensorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_DiagTensor_H
#define Foam_DiagTensor_H

#include "Tensor.H"
#include "SphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DiagTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class DiagTensor
:
    public VectorSpace<DiagTensor<Cmpt>, Cmpt, 3>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef DiagTensor<label> labelType;


    // Member Constants

        //- Rank of DiagTensor is 2
        static constexpr direction rank = 2;


    //- Component labeling enumeration
    enum components { XX, YY, ZZ };


    // Generated Methods

        //- Default construct
        DiagTensor() = default;

        //- Copy construct
        DiagTensor(const DiagTensor&) = default;

        //- Copy assignment
        DiagTensor& operator=(const DiagTensor&) = default;


    // Constructors

        //- Construct initialized to zero
        inline DiagTensor(const Foam::zero);

        //- Construct given VectorSpace
        template<class Cmpt2>
        inline DiagTensor(const VectorSpace<DiagTensor<Cmpt2>, Cmpt2, 3>&);

        //- Construct given SphericalTensor
        inline DiagTensor(const SphericalTensor<Cmpt>&);

        //- Construct given three components
        inline DiagTensor(const Cmpt& txx, const Cmpt& tyy, const Cmpt& tzz);

        //- Construct from Istream
        inline explicit DiagTensor(Istream& is);


    // Member Functions

    // Component Access

        const Cmpt& xx() const noexcept { return this->v_[XX]; }
        const Cmpt& yy() const noexcept { return this->v_[YY]; }
        const Cmpt& zz() const noexcept { return this->v_[ZZ]; }

        Cmpt& xx() noexcept { return this->v_[XX]; }
        Cmpt& yy() noexcept { return this->v_[YY]; }
        Cmpt& zz() noexcept { return this->v_[ZZ]; }


    // Diagonal access and manipulation

        //- The L2-norm squared of the diagonal
        inline scalar diagSqr() const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<DiagTensor<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<DiagTensor<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<DiagTensor<Cmpt>> : is_contiguous_scalar<Cmpt> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DiagTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

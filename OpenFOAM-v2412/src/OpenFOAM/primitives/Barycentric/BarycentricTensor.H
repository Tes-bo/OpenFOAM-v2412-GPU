/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::BarycentricTensor

Description
    Templated 4x3 tensor derived from VectorSpace. Has 12 components.
    Can represent a barycentric transformation as a matrix-barycentric
    inner-product.

SourceFiles
    BarycentricTensorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_BarycentricTensor_H
#define Foam_BarycentricTensor_H

#include "Barycentric.H"
#include "Tensor.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class BarycentricTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class BarycentricTensor
:
    public MatrixSpace<BarycentricTensor<Cmpt>, Cmpt, 4, 3>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        //- (unused)
        typedef BarycentricTensor<label> labelType;


    // Member Constants

        //- Rank of BarycentricTensor is 2
        static constexpr direction rank = 2;


    //- Component labeling enumeration
    enum components { XA, XB, XC, XD, YA, YB, YC, YD, ZA, ZB, ZC, ZD };


    // Generated Methods: copy construct/assignment

        //- Default construct
        BarycentricTensor() = default;


    // Constructors

        //- Construct initialised to zero
        inline BarycentricTensor(const Foam::zero);

        //- Construct given three barycentric components (rows)
        inline BarycentricTensor
        (
            const Barycentric<Cmpt>& x,
            const Barycentric<Cmpt>& y,
            const Barycentric<Cmpt>& z
        );

        //- Construct given four vector components (columns)
        //  Eg, the corners of a tetrahedron
        inline BarycentricTensor
        (
            const Vector<Cmpt>& a,
            const Vector<Cmpt>& b,
            const Vector<Cmpt>& c,
            const Vector<Cmpt>& d
        );


    // Member Functions

    // Component Access

        const Cmpt& xa() const noexcept { return this->v_[XA]; }
        const Cmpt& xb() const noexcept { return this->v_[XB]; }
        const Cmpt& xc() const noexcept { return this->v_[XC]; }
        const Cmpt& xd() const noexcept { return this->v_[XD]; }

        const Cmpt& ya() const noexcept { return this->v_[YA]; }
        const Cmpt& yb() const noexcept { return this->v_[YB]; }
        const Cmpt& yc() const noexcept { return this->v_[YC]; }
        const Cmpt& yd() const noexcept { return this->v_[YD]; }

        const Cmpt& za() const noexcept { return this->v_[ZA]; }
        const Cmpt& zb() const noexcept { return this->v_[ZB]; }
        const Cmpt& zc() const noexcept { return this->v_[ZC]; }
        const Cmpt& zd() const noexcept { return this->v_[ZD]; }


    // Row-barycentric access

        inline Barycentric<Cmpt> x() const;
        inline Barycentric<Cmpt> y() const;
        inline Barycentric<Cmpt> z() const;


    // Column-vector access

        inline Vector<Cmpt> a() const;
        inline Vector<Cmpt> b() const;
        inline Vector<Cmpt> c() const;
        inline Vector<Cmpt> d() const;


    // Operations

        //- Tensor/barycentric inner product
        //  (transforms barycentric coordinates to vector)
        inline Vector<Cmpt> inner(const Barycentric<Cmpt>& bry) const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<BarycentricTensor<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<BarycentricTensor<Cmpt>>
:
    is_contiguous_label<Cmpt>
{};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<BarycentricTensor<Cmpt>>
:
    is_contiguous_scalar<Cmpt>
{};


template<class Cmpt>
class typeOfTranspose<Cmpt, BarycentricTensor<Cmpt>>
{
public:

    typedef void type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BarycentricTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

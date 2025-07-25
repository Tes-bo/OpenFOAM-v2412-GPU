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
    Foam::Barycentric2D

Description
    Templated 2D Barycentric derived from VectorSpace.
    Has 3 components, one of which is redundant.

SourceFiles
    Barycentric2DI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Barycentric2D_H
#define Foam_Barycentric2D_H

#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Barycentric2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Barycentric2D
:
    public VectorSpace<Barycentric2D<Cmpt>, Cmpt, 3>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Barycentric2D<label> labelType;


    // Member Constants

        //- Rank of Barycentric2D is 1
        static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { A, B, C };


    // Generated Methods: copy construct/assignment

        //- Default construct
        Barycentric2D() = default;


    // Constructors

        //- Construct initialized to zero
        inline Barycentric2D(const Foam::zero);

        //- Construct from components
        inline Barycentric2D(const Cmpt& va, const Cmpt& vb, const Cmpt& vc);

        //- Construct from two components, calculating third component
        inline Barycentric2D(const Cmpt& va, const Cmpt& vb);


    // Member Functions

    // Component Access

        const Cmpt& a() const noexcept { return this->v_[A]; }
        const Cmpt& b() const noexcept { return this->v_[B]; }
        const Cmpt& c() const noexcept { return this->v_[C]; }

        Cmpt& a() noexcept { return this->v_[A]; }
        Cmpt& b() noexcept { return this->v_[B]; }
        Cmpt& c() noexcept { return this->v_[C]; }


    // Operations, Tests

        //- Scalar-product of \c this with another Barycentric2D.
        inline Cmpt inner(const Barycentric2D<Cmpt>& b2) const;

        //- True if any coordinate is negative
        inline bool outside() const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Barycentric2D<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Barycentric2D<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Barycentric2D<Cmpt>> : is_contiguous_scalar<Cmpt>{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Barycentric2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::Barycentric

Description
    Templated 3D Barycentric derived from VectorSpace.
    Has 4 components, one of which is redundant.

SourceFiles
    BarycentricI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Barycentric_H
#define Foam_Barycentric_H

#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Barycentric Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Barycentric
:
    public VectorSpace<Barycentric<Cmpt>, Cmpt, 4>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Barycentric<label> labelType;


    // Member Constants

        //- Rank of Barycentric is 1
        static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { A, B, C, D };


    // Generated Methods: copy construct/assignment

        //- Default construct
        Barycentric() = default;


    // Constructors

        //- Construct initialized to zero
        inline Barycentric(const Foam::zero);

        //- Construct from components
        inline Barycentric
        (
            const Cmpt& va,
            const Cmpt& vb,
            const Cmpt& vc,
            const Cmpt& vd
        );

        //- Construct from three components, calculate fourth component
        inline Barycentric(const Cmpt& va, const Cmpt& vb, const Cmpt& vc);


    // Member Functions

    // Component Access

        const Cmpt& a() const noexcept { return this->v_[A]; }
        const Cmpt& b() const noexcept { return this->v_[B]; }
        const Cmpt& c() const noexcept { return this->v_[C]; }
        const Cmpt& d() const noexcept { return this->v_[D]; }

        Cmpt& a() noexcept { return this->v_[A]; }
        Cmpt& b() noexcept { return this->v_[B]; }
        Cmpt& c() noexcept { return this->v_[C]; }
        Cmpt& d() noexcept { return this->v_[D]; }


    // Operations

        //- Scalar-product of \c this with another Barycentric.
        inline Cmpt inner(const Barycentric<Cmpt>& b2) const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Barycentric<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Barycentric<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Barycentric<Cmpt>> : is_contiguous_scalar<Cmpt> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BarycentricI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

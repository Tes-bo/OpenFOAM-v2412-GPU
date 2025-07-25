/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::SpatialVector

Description
    Templated 3D spatial vector derived from VectorSpace used to represent the
    anglular and linear components of position, velocity and acceleration of
    rigid bodies.

    Reference:
    \verbatim
        Featherstone, R. (2008).
        Rigid body dynamics algorithms.
        Springer.
    \endverbatim

SourceFiles
    SpatialVectorI.H

See also
    Foam::VectorSpace
    Foam::Vector

\*---------------------------------------------------------------------------*/

#ifndef Foam_SpatialVector_H
#define Foam_SpatialVector_H

#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SpatialVector Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SpatialVector
:
    public VectorSpace<SpatialVector<Cmpt>, Cmpt, 6>
{

public:

    //- Component labeling enumeration
    enum components { WX, WY, WZ, LX, LY, LZ };


    //- Class to represent the dual spatial vector
    class dual
    {
        const SpatialVector& v_;

    public:

        //- Construct the dual of the given SpatialVector
        inline dual(const SpatialVector& v);

        //- Return the parent SpatialVector
        inline const SpatialVector& v() const;
    };


    // Constructors

        //- Default construct
        SpatialVector() = default;

        //- Construct initialized to zero
        inline SpatialVector(const Foam::zero);

        //- Construct given VectorSpace of the same rank
        inline SpatialVector(const typename SpatialVector::vsType&);

        //- Construct from the angular and linear vector components
        inline SpatialVector
        (
            const Vector<Cmpt>& w,
            const Vector<Cmpt>& l
        );

        //- Construct given 6 components
        inline SpatialVector
        (
            const Cmpt& wx,
            const Cmpt& wy,
            const Cmpt& wz,
            const Cmpt& lx,
            const Cmpt& ly,
            const Cmpt& lz
        );

        //- Construct from Istream
        inline explicit SpatialVector(Istream&);


    // Member Functions

    // Component Access

        const Cmpt& wx() const noexcept { return this->v_[WX]; }
        const Cmpt& wy() const noexcept { return this->v_[WY]; }
        const Cmpt& wz() const noexcept { return this->v_[WZ]; }

        const Cmpt& lx() const noexcept { return this->v_[LX]; }
        const Cmpt& ly() const noexcept { return this->v_[LY]; }
        const Cmpt& lz() const noexcept { return this->v_[LZ]; }

        Cmpt& wx() noexcept { return this->v_[WX]; }
        Cmpt& wy() noexcept { return this->v_[WY]; }
        Cmpt& wz() noexcept { return this->v_[WZ]; }

        Cmpt& lx() noexcept { return this->v_[LX]; }
        Cmpt& ly() noexcept { return this->v_[LY]; }
        Cmpt& lz() noexcept { return this->v_[LZ]; }


    // Sub-vector access

        //- Return the angular part of the spatial vector as a vector
        inline Vector<Cmpt> w() const;

        //- Return the linear part of the spatial vector as a vector
        inline Vector<Cmpt> l() const;


    // Member Operators

        //- Return the dual spatial vector
        inline dual operator*() const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<SpatialVector<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<SpatialVector<Cmpt>>
:
    is_contiguous_label<Cmpt>
{};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<SpatialVector<Cmpt>>
:
    is_contiguous_scalar<Cmpt>
{};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "SpatialVectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Vector<Cmpt>::Vector(const Foam::zero)
:
    Vector::vsType(Zero)
{}


template<class Cmpt>
template<class Cmpt2>
inline Foam::Vector<Cmpt>::Vector
(
    const VectorSpace<Vector<Cmpt2>, Cmpt2, 3>& vs
)
:
    Vector::vsType(vs)
{}


template<class Cmpt>
inline Foam::Vector<Cmpt>::Vector
(
    const Cmpt& vx,
    const Cmpt& vy,
    const Cmpt& vz
)
{
    this->v_[X] = vx;
    this->v_[Y] = vy;
    this->v_[Z] = vz;
}


template<class Cmpt>
inline Foam::Vector<Cmpt>::Vector(Istream& is)
:
    Vector::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Foam::Vector<Cmpt>& Foam::Vector<Cmpt>::centre
(
    const Foam::UList<Vector<Cmpt>>&
) const noexcept
{
    return *this;
}


template<class Cmpt>
inline Foam::scalar Foam::Vector<Cmpt>::magSqr() const
{
    return
    (
        Foam::magSqr(this->x())
      + Foam::magSqr(this->y())
      + Foam::magSqr(this->z())
    );
}


template<class Cmpt>
inline Foam::scalar Foam::Vector<Cmpt>::mag() const
{
    return ::sqrt(this->magSqr());
}


template<class Cmpt>
inline Foam::scalar Foam::Vector<Cmpt>::distSqr(const Vector<Cmpt>& v2) const
{
    return
    (
        Foam::magSqr(v2.x() - this->x())
      + Foam::magSqr(v2.y() - this->y())
      + Foam::magSqr(v2.z() - this->z())
    );
}


template<class Cmpt>
inline Foam::scalar Foam::Vector<Cmpt>::dist(const Vector<Cmpt>& v2) const
{
    return ::sqrt(this->distSqr(v2));
}


template<class Cmpt>
inline Foam::Vector<Cmpt>& Foam::Vector<Cmpt>::normalise(const scalar tol)
{
    #ifdef __clang__
    volatile  // Use volatile to avoid aggressive branch optimization
    #endif
    const scalar s = this->mag();

    if (s < tol)
    {
        *this = Zero;
    }
    else
    {
        *this /= s;
    }

    return *this;
}


template<class Cmpt>
inline Foam::Vector<Cmpt>&
Foam::Vector<Cmpt>::removeCollinear(const Vector<Cmpt>& unitVec)
{
    *this -= (*this & unitVec) * unitVec;
    return *this;
}


// * * * * * * * * * * * * * * * Member Operations * * * * * * * * * * * * * //

template<class Cmpt>
inline Cmpt Foam::Vector<Cmpt>::inner(const Vector<Cmpt>& v2) const
{
    const Vector<Cmpt>& v1 = *this;

    return (v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z());
}


template<class Cmpt>
inline Foam::Vector<Cmpt>
Foam::Vector<Cmpt>::cross(const Vector<Cmpt>& v2) const
{
    const Vector<Cmpt>& v1 = *this;

    return Vector<Cmpt>
    (
        (v1.y()*v2.z() - v1.z()*v2.y()),
        (v1.z()*v2.x() - v1.x()*v2.z()),
        (v1.x()*v2.y() - v1.y()*v2.x())
    );
}


// * * * * * * * * * * * * * Comparison Operations * * * * * * * * * * * * * //

template<class Cmpt>
inline bool
Foam::Vector<Cmpt>::less_xyz(const Vector<Cmpt>& a, const Vector<Cmpt>& b)
{
    return
    (
        (a.x() < b.x())  // Component is less
     ||
        (
            !(b.x() < a.x())  // Equal? Check next component
         &&
            (
                (a.y() < b.y())  // Component is less
             ||
                (
                    !(b.y() < a.y())  // Equal? Check next component
                 && (a.z() < b.z())
                )
            )
        )
    );
}


template<class Cmpt>
inline bool
Foam::Vector<Cmpt>::less_yzx(const Vector<Cmpt>& a, const Vector<Cmpt>& b)
{
    return
    (
        (a.y() < b.y())  // Component is less
     ||
        (
            !(b.y() < a.y())  // Equal? Check next component
         &&
            (
                (a.z() < b.z())  // Component is less
             ||
                (
                    !(b.z() < a.z())  // Equal? Check next component
                 && (a.x() < b.x())
                )
            )
        )
    );
}


template<class Cmpt>
inline bool
Foam::Vector<Cmpt>::less_zxy(const Vector<Cmpt>& a, const Vector<Cmpt>& b)
{
    return
    (
        (a.z() < b.z())  // Component is less
     ||
        (
            !(b.z() < a.z())  // Equal? Check next component
         &&
            (
                (a.x() < b.x())  // Component is less
             ||
                (
                    !(b.x() < a.x())  // Equal? Check next component
                 && (a.y() < b.y())
                )
            )
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Linear interpolation of vectors a and b by factor t
template<class Cmpt>
inline Vector<Cmpt> lerp
(
    const Vector<Cmpt>& a,
    const Vector<Cmpt>& b,
    const scalar t
)
{
    const scalar onet = (1-t);

    return Vector<Cmpt>
    (
        onet*a.x() + t*b.x(),
        onet*a.y() + t*b.y(),
        onet*a.z() + t*b.z()
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Dummy innerProduct for scalar
//  Allows the construction of vtables for virtual member functions
//  involving the inner-products of fields
//  for which a "NotImplemented" specialization for scalar is provided.
template<class Cmpt>
class innerProduct<Vector<Cmpt>, scalar>
{
public:

    typedef scalar type;
};


template<class Cmpt>
inline Cmpt operator&(const Vector<Cmpt>& v1, const Vector<Cmpt>& v2)
{
    return v1.inner(v2);
}


template<class Cmpt>
inline Vector<Cmpt> operator^(const Vector<Cmpt>& v1, const Vector<Cmpt>& v2)
{
    return v1.cross(v2);
}


template<class Cmpt>
inline Vector<Cmpt> operator*(const Cmpt& s, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>(s*v.x(), s*v.y(), s*v.z());
}


template<class Cmpt>
inline Vector<Cmpt> operator*(const Vector<Cmpt>& v, const Cmpt& s)
{
    return s*v;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

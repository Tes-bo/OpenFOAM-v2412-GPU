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
inline Foam::Tensor2D<Cmpt>::Tensor2D(const Foam::zero)
:
    Tensor2D::vsType(Zero)
{}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>::Tensor2D
(
    const VectorSpace<Tensor2D<Cmpt>, Cmpt, 4>& vs
)
:
    Tensor2D::vsType(vs)
{}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>::Tensor2D(const SymmTensor2D<Cmpt>& st)
{
    this->v_[XX] = st.xx(); this->v_[XY] = st.xy();
    this->v_[YX] = st.xy(); this->v_[YY] = st.yy();
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>::Tensor2D(const SphericalTensor2D<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = Zero;
    this->v_[YX] = Zero;    this->v_[YY] = st.ii();
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>::Tensor2D
(
    const Vector2D<Cmpt>& x,
    const Vector2D<Cmpt>& y
)
{
    this->v_[XX] = x.x(); this->v_[XY] = x.y();
    this->v_[YX] = y.x(); this->v_[YY] = y.y();
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>::Tensor2D
(
    const Cmpt txx, const Cmpt txy,
    const Cmpt tyx, const Cmpt tyy
)
{
    this->v_[XX] = txx; this->v_[XY] = txy;
    this->v_[YX] = tyx; this->v_[YY] = tyy;
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>::Tensor2D(Istream& is)
:
    Tensor2D::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::x() const
{
    return Vector2D<Cmpt>(this->v_[XX], this->v_[XY]);
}

template<class Cmpt>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::y() const
{
    return Vector2D<Cmpt>(this->v_[YX], this->v_[YY]);
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::cx() const
{
    return Vector2D<Cmpt>(this->v_[XX], this->v_[YX]);
}

template<class Cmpt>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::cy() const
{
    return Vector2D<Cmpt>(this->v_[XY], this->v_[YY]);
}


template<class Cmpt>
template<Foam::direction Idx>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::col() const
{
    if (Idx == 0) return cx();
    else if (Idx == 1) return cy();

    static_assert(Idx < 2, "Invalid column access");
    return Zero;
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::col(const direction c) const
{
    switch (c)
    {
        case 0: return cx(); break;
        case 1: return cy(); break;
        default:
            FatalErrorInFunction
                << "Invalid column access " << c << abort(FatalError);
    }

    return Zero;
}


template<class Cmpt>
template<Foam::direction Idx>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::row() const
{
    if (Idx == 0) return x();
    else if (Idx == 1) return y();

    static_assert(Idx < 2, "Invalid row access");
    return Zero;
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::row(const direction r) const
{
    switch (r)
    {
        case 0: return x(); break;
        case 1: return y(); break;
        default:
            FatalErrorInFunction
                << "Invalid row access " << r << abort(FatalError);
    }

    return Zero;
}


template<class Cmpt>
template<Foam::direction Idx>
inline void Foam::Tensor2D<Cmpt>::col(const Vector2D<Cmpt>& v)
{
    if (Idx == 0)
    {
        this->v_[XX] = v.x();
        this->v_[YX] = v.y();
    }
    else if (Idx == 1)
    {
        this->v_[XY] = v.x();
        this->v_[YY] = v.y();
    }

    static_assert(Idx < 2, "Invalid column access");
}


template<class Cmpt>
template<Foam::direction Idx>
inline void Foam::Tensor2D<Cmpt>::row(const Vector2D<Cmpt>& v)
{
    if (Idx == 0)
    {
        this->v_[XX] = v.x(); this->v_[XY] = v.y();
    }
    else if (Idx == 1)
    {
        this->v_[YX] = v.x(); this->v_[YY] = v.y();
    }

    static_assert(Idx < 2, "Invalid row access");
}


template<class Cmpt>
inline void Foam::Tensor2D<Cmpt>::cols
(
    const Vector2D<Cmpt>& x,
    const Vector2D<Cmpt>& y
)
{
    this->v_[XX] = x.x(); this->v_[XY] = y.x();
    this->v_[YX] = x.y(); this->v_[YY] = y.y();
}


template<class Cmpt>
inline void Foam::Tensor2D<Cmpt>::rows
(
    const Vector2D<Cmpt>& x,
    const Vector2D<Cmpt>& y
)
{
    this->v_[XX] = x.x(); this->v_[XY] = x.y();
    this->v_[YX] = y.x(); this->v_[YY] = y.y();
}


template<class Cmpt>
inline void Foam::Tensor2D<Cmpt>::col
(
    const direction c,
    const Vector2D<Cmpt>& v
)
{
    switch (c)
    {
        case 0: col<0>(v); break;
        case 1: col<1>(v); break;
        default:
            FatalErrorInFunction
                << "Invalid column access " << c << abort(FatalError);
    }
}


template<class Cmpt>
inline void Foam::Tensor2D<Cmpt>::row
(
    const direction r,
    const Vector2D<Cmpt>& v
)
{
    switch (r)
    {
        case 0: row<0>(v); break;
        case 1: row<1>(v); break;
        default:
            FatalErrorInFunction
                << "Invalid row access " << r << abort(FatalError);
    }
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt> Foam::Tensor2D<Cmpt>::diag() const
{
    return Vector2D<Cmpt>(this->v_[XX], this->v_[YY]);
}


template<class Cmpt>
inline void Foam::Tensor2D<Cmpt>::diag(const Vector2D<Cmpt>& v)
{
    this->v_[XX] = v.x(); this->v_[YY] = v.y();
}


template<class Cmpt>
inline Foam::scalar Foam::Tensor2D<Cmpt>::diagSqr() const
{
    return
    (
        Foam::magSqr(this->xx())
      + Foam::magSqr(this->yy())
    );
}


// * * * * * * * * * * * * * * * Member Operations * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Tensor2D<Cmpt> Foam::Tensor2D<Cmpt>::T() const
{
    return Tensor2D<Cmpt>
    (
        xx(), yx(),
        xy(), yy()
    );
}


template<class Cmpt>
inline Cmpt Foam::Tensor2D<Cmpt>::det() const
{
    return (xx()*yy() - xy()*yx());
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt> Foam::Tensor2D<Cmpt>::adjunct() const
{
    return Tensor2D<Cmpt>
    (
        yy(), -xy(),
       -yx(),  xx()
    );
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt> Foam::Tensor2D<Cmpt>::cof() const
{
    return this->adjunct().T();
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>
Foam::Tensor2D<Cmpt>::inner(const Tensor2D<Cmpt>& t2) const
{
    const Tensor2D<Cmpt>& t1 = *this;

    return Tensor2D<Cmpt>
    (
        t1.xx()*t2.xx() + t1.xy()*t2.yx(),
        t1.xx()*t2.xy() + t1.xy()*t2.yy(),

        t1.yx()*t2.xx() + t1.yy()*t2.yx(),
        t1.yx()*t2.xy() + t1.yy()*t2.yy()
    );
}


template<class Cmpt>
inline Foam::Tensor2D<Cmpt>
Foam::Tensor2D<Cmpt>::schur(const Tensor2D<Cmpt>& t2) const
{
    const Tensor2D<Cmpt>& t1 = *this;

    return Tensor2D<Cmpt>
    (
        t1.xx()*t2.xx(), t1.xy()*t2.xy(),
        t1.yx()*t2.yx(), t1.yy()*t2.yy()
    );
}


// Invert without much error handling
template<class Cmpt>
inline Foam::Tensor2D<Cmpt> Foam::Tensor2D<Cmpt>::inv() const
{
    const Cmpt detval = this->det();

    #ifdef FULLDEBUG
    if (mag(detval) < SMALL)
    {
        FatalErrorInFunction
            << "SymmTensor2D not properly invertible, determinant:"
            << detval << " tensor:" << *this << nl
            << abort(FatalError);
    }
    #endif

    return this->adjunct()/detval;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline void Foam::Tensor2D<Cmpt>::operator=(const SymmTensor2D<Cmpt>& st)
{
    this->v_[XX] = st.xx(); this->v_[XY] = st.xy();
    this->v_[YX] = st.xy(); this->v_[YY] = st.yy();
}


template<class Cmpt>
inline void Foam::Tensor2D<Cmpt>::operator=(const SphericalTensor2D<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = Zero;
    this->v_[YX] = Zero;    this->v_[YY] = st.ii();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the trace of a Tensor2D
template<class Cmpt>
inline Cmpt tr(const Tensor2D<Cmpt>& t)
{
    return t.xx() + t.yy();
}


//- Return the spherical part of a Tensor2D
template<class Cmpt>
inline SphericalTensor2D<Cmpt> sph(const Tensor2D<Cmpt>& t)
{
    return SphericalTensor2D<Cmpt>
    (
        0.5*tr(t)
    );
}


//- Return the symmetric part of a Tensor2D
template<class Cmpt>
inline SymmTensor2D<Cmpt> symm(const Tensor2D<Cmpt>& t)
{
    return SymmTensor2D<Cmpt>
    (
        t.xx(), 0.5*(t.xy() + t.yx()),
                t.yy()
    );
}


//- Return the twice the symmetric part of a Tensor2D
template<class Cmpt>
inline SymmTensor2D<Cmpt> twoSymm(const Tensor2D<Cmpt>& t)
{
    return SymmTensor2D<Cmpt>
    (
        t.xx() + t.xx(), t.xy() + t.yx(),
                         t.yy() + t.yy()
    );
}


//- Return the skew-symmetric part of a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt> skew(const Tensor2D<Cmpt>& t)
{
    return Tensor2D<Cmpt>
    (
        Zero,                  0.5*(t.xy() - t.yx()),
        0.5*(t.yx() - t.xy()), Zero
    );
}


//- Return the deviatoric part of a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt> dev(const Tensor2D<Cmpt>& t)
{
    return t - sph(t);
}


//- Return the two-third deviatoric part of a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt> dev2(const Tensor2D<Cmpt>& t)
{
    return t - 2*sph(t);
}


//- Return the determinant of a Tensor2D
template<class Cmpt>
inline Cmpt det(const Tensor2D<Cmpt>& t)
{
    return t.det();
}


//- Return the cofactor of a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt> cof(const Tensor2D<Cmpt>& t)
{
    return t.cof();
}


//- Return the inverse of a Tensor2D, using given determinant value
template<class Cmpt>
inline Tensor2D<Cmpt> inv(const Tensor2D<Cmpt>& t, const Cmpt detval)
{
    #ifdef FULLDEBUG
    if (mag(detval) < SMALL)
    {
        FatalErrorInFunction
            << "Tensor2D not properly invertible, determinant:"
            << detval << " tensor:" << t << nl
            << abort(FatalError);
    }
    #endif

    return t.adjunct()/detval;
}


//- Return the inverse of a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt> inv(const Tensor2D<Cmpt>& t)
{
    return t.inv();
}


//- Return the 1st invariant of a Tensor2D
template<class Cmpt>
inline Cmpt invariantI(const Tensor2D<Cmpt>& t)
{
    return tr(t);
}


//- Return the 2nd invariant of a Tensor2D
template<class Cmpt>
inline Cmpt invariantII(const Tensor2D<Cmpt>& t)
{
    return det(t);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Sum of a SphericalTensor2D and a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator+(const SphericalTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return Tensor2D<Cmpt>
    (
        st1.ii() + t2.xx(), t2.xy(),
        t2.yx(),            st1.ii() + t2.yy()
    );
}


//- Sum of a Tensor2D and a SphericalTensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator+(const Tensor2D<Cmpt>& t1, const SphericalTensor2D<Cmpt>& st2)
{
    return Tensor2D<Cmpt>
    (
        t1.xx() + st2.ii(), t1.xy(),
        t1.yx(),            t1.yy() + st2.ii()
    );
}


//- Sum of a SymmTensor2D and a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator+(const SymmTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return Tensor2D<Cmpt>
    (
        st1.xx() + t2.xx(), st1.xy() + t2.xy(),
        st1.xy() + t2.yx(), st1.yy() + t2.yy()
    );
}


//- Sum of a Tensor2D and a SymmTensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator+(const Tensor2D<Cmpt>& t1, const SymmTensor2D<Cmpt>& st2)
{
    return Tensor2D<Cmpt>
    (
        t1.xx() + st2.xx(), t1.xy() + st2.xy(),
        t1.yx() + st2.xy(), t1.yy() + st2.yy()
    );
}


//- Subtract a Tensor2D from a SphericalTensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator-(const SphericalTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return Tensor2D<Cmpt>
    (
        st1.ii() - t2.xx(), -t2.xy(),
       -t2.yx(),             st1.ii() - t2.yy()
    );
}


//- Subtract a SphericalTensor2D from a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator-(const Tensor2D<Cmpt>& t1, const SphericalTensor2D<Cmpt>& st2)
{
    return Tensor2D<Cmpt>
    (
        t1.xx() - st2.ii(), t1.xy(),
        t1.yx(),            t1.yy() - st2.ii()
    );
}


//- Subtract a Tensor2D from a SymmTensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator-(const SymmTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return Tensor2D<Cmpt>
    (
        st1.xx() - t2.xx(), st1.xy() - t2.xy(),
        st1.xy() - t2.yx(), st1.yy() - t2.yy()
    );
}


//- Subtract a SymmTensor2D from a Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator-(const Tensor2D<Cmpt>& t1, const SymmTensor2D<Cmpt>& st2)
{
    return Tensor2D<Cmpt>
    (
        t1.xx() - st2.xx(), t1.xy() - st2.xy(),
        t1.yx() - st2.xy(), t1.yy() - st2.yy()
    );
}


//- Division of a Tensor2D by a Cmpt
template<class Cmpt>
inline Tensor2D<Cmpt>
operator/(const Tensor2D<Cmpt>& t, const Cmpt s)
{
    #ifdef FULLDEBUG
    if (mag(s) < VSMALL)
    {
        FatalErrorInFunction
            << "Tensor2D = " << t
            << " is not divisible due to a zero value in Cmpt:"
            << "Cmpt = " << s
            << abort(FatalError);
    }
    #endif

    return Tensor2D<Cmpt>
    (
        t.xx()/s, t.xy()/s,
        t.yx()/s, t.yy()/s
    );
}


//- Inner-product of a Tensor2D and a Tensor2D
template<class Cmpt>
inline typename innerProduct<Tensor2D<Cmpt>, Tensor2D<Cmpt>>::type
operator&(const Tensor2D<Cmpt>& t1, const Tensor2D<Cmpt>& t2)
{
    return t1.inner(t2);
}


//- Inner-product of a SphericalTensor2D and Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator&(const SphericalTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return Tensor2D<Cmpt>
    (
        st1.ii()*t2.xx(),
        st1.ii()*t2.xy(),
                          st1.ii()*t2.yx(),
                          st1.ii()*t2.yy()
    );
}


//- Inner-product of a Tensor2D and a SphericalTensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator&(const Tensor2D<Cmpt>& t1, const SphericalTensor2D<Cmpt>& st2)
{
    return Tensor2D<Cmpt>
    (
        t1.xx()*st2.ii(),
                          t1.xy()*st2.ii(),

        t1.yx()*st2.ii(),
                          t1.yy()*st2.ii()
    );
}


//- Inner-product of a SymmTensor2D and Tensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator&(const SymmTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return Tensor2D<Cmpt>
    (
        st1.xx()*t2.xx() + st1.xy()*t2.yx(),
        st1.xx()*t2.xy() + st1.xy()*t2.yy(),

        st1.xy()*t2.xx() + st1.yy()*t2.yx(),
        st1.xy()*t2.xy() + st1.yy()*t2.yy()
    );
}


//- Inner-product of a Tensor2D and a SymmTensor2D
template<class Cmpt>
inline Tensor2D<Cmpt>
operator&(const Tensor2D<Cmpt>& t1, const SymmTensor2D<Cmpt>& st2)
{
    return Tensor2D<Cmpt>
    (
        t1.xx()*st2.xx() + t1.xy()*st2.xy(),
        t1.xx()*st2.xy() + t1.xy()*st2.yy(),

        t1.yx()*st2.xx() + t1.yy()*st2.xy(),
        t1.yx()*st2.xy() + t1.yy()*st2.yy()
    );
}



//- Inner-product of a Tensor2D and a Vector2D
template<class Cmpt>
inline typename innerProduct<Tensor2D<Cmpt>, Vector2D<Cmpt>>::type
operator&(const Tensor2D<Cmpt>& t, const Vector2D<Cmpt>& v)
{
    return Vector2D<Cmpt>
    (
        t.xx()*v.x() + t.xy()*v.y(),
        t.yx()*v.x() + t.yy()*v.y()
    );
}


//- Inner-product of a Vector2D and a Tensor2D
template<class Cmpt>
inline typename innerProduct<Vector2D<Cmpt>, Tensor2D<Cmpt>>::type
operator&(const Vector2D<Cmpt>& v, const Tensor2D<Cmpt>& t)
{
    return Vector2D<Cmpt>
    (
        v.x()*t.xx() + v.y()*t.yx(),
        v.x()*t.xy() + v.y()*t.yy()
    );
}


//- Double-inner-product of a SphericalTensor2D and a Tensor2D
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return (st1.ii()*t2.xx() + st1.ii()*t2.yy());
}


//- Double-inner-product of a Tensor2D and a SphericalTensor2D
template<class Cmpt>
inline Cmpt
operator&&(const Tensor2D<Cmpt>& t1, const SphericalTensor2D<Cmpt>& st2)
{
    return (t1.xx()*st2.ii() + t1.yy()*st2.ii());
}


//- Double-inner-product of a SymmTensor2D and a Tensor2D
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor2D<Cmpt>& st1, const Tensor2D<Cmpt>& t2)
{
    return
    (
        st1.xx()*t2.xx() + st1.xy()*t2.xy()
      + st1.xy()*t2.yx() + st1.yy()*t2.yy()
    );
}


//- Double-inner-product of a Tensor2D and a SymmTensor2D
template<class Cmpt>
inline Cmpt
operator&&(const Tensor2D<Cmpt>& t1, const SymmTensor2D<Cmpt>& st2)
{
    return
    (
        t1.xx()*st2.xx() + t1.xy()*st2.xy()
      + t1.yx()*st2.xy() + t1.yy()*st2.yy()
    );
}


//- Outer-product of a Vector2D and a Vector2D
template<class Cmpt>
inline typename outerProduct<Vector2D<Cmpt>, Vector2D<Cmpt>>::type
operator*(const Vector2D<Cmpt>& v1, const Vector2D<Cmpt>& v2)
{
    return Tensor2D<Cmpt>
    (
        v1.x()*v2.x(), v1.x()*v2.y(),
        v1.y()*v2.x(), v1.y()*v2.y()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class typeOfSum<SphericalTensor2D<Cmpt>, Tensor2D<Cmpt>>
{
public:

    typedef Tensor2D<Cmpt> type;
};

template<class Cmpt>
class typeOfSum<Tensor2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef Tensor2D<Cmpt> type;
};


template<class Cmpt>
class innerProduct<Tensor2D<Cmpt>, Tensor2D<Cmpt>>
{
public:

    typedef Tensor2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SphericalTensor2D<Cmpt>, Tensor2D<Cmpt>>
{
public:

    typedef Tensor2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Tensor2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef Tensor2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Tensor2D<Cmpt>, Vector2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector2D<Cmpt>, Tensor2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};


template<class Cmpt>
class outerProduct<Vector2D<Cmpt>, Vector2D<Cmpt>>
{
public:

    typedef Tensor2D<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "Vector.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::SphericalTensor<Cmpt>::SphericalTensor(const Foam::zero)
:
    SphericalTensor::vsType(Zero)
{}


template<class Cmpt>
template<class Cmpt2>
inline Foam::SphericalTensor<Cmpt>::SphericalTensor
(
    const VectorSpace<SphericalTensor<Cmpt2>, Cmpt2, 1>& vs
)
:
    SphericalTensor::vsType(vs)
{}


template<class Cmpt>
inline Foam::SphericalTensor<Cmpt>::SphericalTensor(const Cmpt& stii)
{
    this->v_[II] = stii;
}


template<class Cmpt>
inline Foam::SphericalTensor<Cmpt>::SphericalTensor(Istream& is)
:
    SphericalTensor::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::scalar Foam::SphericalTensor<Cmpt>::diagSqr() const
{
    return 3*Foam::magSqr(this->ii());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the trace of a SphericalTensor
template<class Cmpt>
inline Cmpt tr(const SphericalTensor<Cmpt>& st)
{
    return 3*st.ii();
}


//- Return the spherical part of a SphericalTensor, i.e. itself
template<class Cmpt>
inline SphericalTensor<Cmpt> sph(const SphericalTensor<Cmpt>& st)
{
    return st;
}


//- Return the determinant of a SphericalTensor
template<class Cmpt>
inline Cmpt det(const SphericalTensor<Cmpt>& st)
{
    return st.ii()*st.ii()*st.ii();
}


//- Return the inverse of a SphericalTensor
template<class Cmpt>
inline SphericalTensor<Cmpt> inv(const SphericalTensor<Cmpt>& st)
{
    #ifdef FULLDEBUG
    if (mag(st.ii()) < VSMALL)
    {
        FatalErrorInFunction
            << "SphericalTensor not invertible, determinant:"
            << det(st) << " tensor:" << st << nl
            << abort(FatalError);
    }
    #endif

    return SphericalTensor<Cmpt>(1/st.ii());
}


//- Return the square of Frobenius norm of a SphericalTensor
template<class Cmpt>
inline Foam::scalar magSqr(const SphericalTensor<Cmpt>& st)
{
    return st.diagSqr();
}


//- Return the max component of a SphericalTensor
template<class Cmpt>
inline Cmpt cmptMax(const SphericalTensor<Cmpt>& st)
{
    return st.ii();
}


//- Return the min component of a SphericalTensor
template<class Cmpt>
inline Cmpt cmptMin(const SphericalTensor<Cmpt>& st)
{
    return st.ii();
}


//- Return the sum of components of a SphericalTensor
template<class Cmpt>
inline Cmpt cmptSum(const SphericalTensor<Cmpt>& st)
{
    return 3*st.ii();
}


//- Return the arithmetic average of components of a SphericalTensor
template<class Cmpt>
inline Cmpt cmptAv(const SphericalTensor<Cmpt>& st)
{
    return st.ii();
}


//- Linear interpolation of spherical tensors a and b by factor t
template<class Cmpt>
inline SphericalTensor<Cmpt> lerp
(
    const SphericalTensor<Cmpt>& a,
    const SphericalTensor<Cmpt>& b,
    const scalar t
)
{
    return SphericalTensor<Cmpt>((1-t)*a.ii() + t*b.ii());
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Division of a Cmpt by a SphericalTensor
template<class Cmpt>
inline SphericalTensor<Cmpt>
operator/(const Cmpt s, const SphericalTensor<Cmpt>& st)
{
    #ifdef FULLDEBUG
    if (mag(st.ii()) < VSMALL)
    {
        FatalErrorInFunction
            << "Cmpt = " << s
            << " is not divisible due to a zero element in SphericalTensor:"
            << "SphericalTensor = " << st
            << abort(FatalError);
    }
    #endif

    return SphericalTensor<Cmpt>(s/st.ii());
}


//- Division of a SphericalTensor by a Cmpt
template<class Cmpt>
inline SphericalTensor<Cmpt>
operator/(const SphericalTensor<Cmpt>& st, const Cmpt s)
{
    #ifdef FULLDEBUG
    if (mag(s) < VSMALL)
    {
        FatalErrorInFunction
            << "SphericalTensor = " << st
            << " is not divisible due to a zero value in Cmpt:"
            << "Cmpt = " << s
            << abort(FatalError);
    }
    #endif

    return SphericalTensor<Cmpt>(st.ii()/s);
}


//- Inner-product of a SphericalTensor and a SphericalTensor
template<class Cmpt>
inline SphericalTensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& st2)
{
    return SphericalTensor<Cmpt>(st1.ii()*st2.ii());
}


//- Inner-product of a SphericalTensor and a Vector
template<class Cmpt>
inline Vector<Cmpt>
operator&(const SphericalTensor<Cmpt>& st, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        st.ii()*v.x(),
                       st.ii()*v.y(),
                                      st.ii()*v.z()
    );
}


//- Inner-product of a Vector and a SphericalTensor
template<class Cmpt>
inline Vector<Cmpt>
operator&(const Vector<Cmpt>& v, const SphericalTensor<Cmpt>& st)
{
    return Vector<Cmpt>
    (
        v.x()*st.ii(),
                       v.y()*st.ii(),
                                      v.z()*st.ii()
    );
}


//- Double-inner-product of a SphericalTensor and a SphericalTensor
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& st2)
{
    return 3*st1.ii()*st2.ii();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class outerProduct<SphericalTensor<Cmpt>, Cmpt>
{
public:

    typedef SphericalTensor<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SphericalTensor<Cmpt>>
{
public:

    typedef SphericalTensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef SphericalTensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, Vector<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

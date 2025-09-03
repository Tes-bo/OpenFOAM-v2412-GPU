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

#include "Vector2D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D(const Foam::zero)
:
    SphericalTensor2D::vsType(Zero)
{}


template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D
(
    const VectorSpace<SphericalTensor2D<Cmpt>, Cmpt, 1>& vs
)
:
    SphericalTensor2D::vsType(vs)
{}


template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D(const Cmpt& stii)
{
    this->v_[II] = stii;
}


template<class Cmpt>
inline Foam::SphericalTensor2D<Cmpt>::SphericalTensor2D(Istream& is)
:
    SphericalTensor2D::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::scalar Foam::SphericalTensor2D<Cmpt>::diagSqr() const
{
    return 2*Foam::magSqr(this->ii());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the trace of a SphericalTensor2D
template<class Cmpt>
inline Cmpt tr(const SphericalTensor2D<Cmpt>& st)
{
    return 2*st.ii();
}


//- Return the spherical part of a SphericalTensor2D, i.e. itself
template<class Cmpt>
inline SphericalTensor2D<Cmpt> sph(const SphericalTensor2D<Cmpt>& st)
{
    return st;
}


//- Return the determinant of a SphericalTensor2D
template<class Cmpt>
inline Cmpt det(const SphericalTensor2D<Cmpt>& st)
{
    return st.ii()*st.ii();
}


//- Return the inverse of a SphericalTensor2D
template<class Cmpt>
inline SphericalTensor2D<Cmpt> inv(const SphericalTensor2D<Cmpt>& st)
{
    #ifdef FULLDEBUG
    if (mag(st.ii()) < VSMALL)
    {
        FatalErrorInFunction
            << "SphericalTensor2D not invertible, determinant:"
            << det(st) << " tensor:" << st << nl
            << abort(FatalError);
    }
    #endif

    return SphericalTensor2D<Cmpt>(1/st.ii());
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Division of a Cmpt by a SphericalTensor2D
template<class Cmpt>
inline SphericalTensor2D<Cmpt>
operator/(const Cmpt s, const SphericalTensor2D<Cmpt>& st)
{
    #ifdef FULLDEBUG
    if (mag(st.ii()) < VSMALL)
    {
        FatalErrorInFunction
            << "Cmpt = " << s
            << " is not divisible due to a zero element in SphericalTensor2D:"
            << "SphericalTensor2D = " << st
            << abort(FatalError);
    }
    #endif

    return SphericalTensor2D<Cmpt>(s/st.ii());
}


//- Division of a SphericalTensor2D by a Cmpt
template<class Cmpt>
inline SphericalTensor2D<Cmpt>
operator/(const SphericalTensor2D<Cmpt>& st, const Cmpt s)
{
    #ifdef FULLDEBUG
    if (mag(s) < VSMALL)
    {
        FatalErrorInFunction
            << "SphericalTensor2D = " << st
            << " is not divisible due to a zero value in Cmpt:"
            << "Cmpt = " << s
            << abort(FatalError);
    }
    #endif

    return SphericalTensor2D<Cmpt>(st.ii()/s);
}


//- Inner-product of a SphericalTensor2D and a SphericalTensor2D
template<class Cmpt>
inline SphericalTensor2D<Cmpt>
operator&
(
    const SphericalTensor2D<Cmpt>& st1,
    const SphericalTensor2D<Cmpt>& st2
)
{
    return SphericalTensor2D<Cmpt>(st1.ii()*st2.ii());
}


//- Inner-product of a SphericalTensor2D and a Vector2D
template<class Cmpt>
inline Vector2D<Cmpt>
operator&(const SphericalTensor2D<Cmpt>& st, const Vector2D<Cmpt>& v)
{
    return Vector2D<Cmpt>
    (
        st.ii()*v.x(),
                       st.ii()*v.y()
    );
}


//- Inner-product of a Vector2D and a SphericalTensor2D
template<class Cmpt>
inline Vector2D<Cmpt>
operator&(const Vector2D<Cmpt>& v, const SphericalTensor2D<Cmpt>& st)
{
    return Vector2D<Cmpt>
    (
        v.x()*st.ii(),
                       v.y()*st.ii()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class outerProduct<SphericalTensor2D<Cmpt>, Cmpt>
{
public:

    typedef SphericalTensor2D<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SphericalTensor2D<Cmpt>>
{
public:

    typedef SphericalTensor2D<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef SphericalTensor2D<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor2D<Cmpt>, Vector2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector2D<Cmpt>, SphericalTensor2D<Cmpt>>
{
public:

    typedef Vector2D<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

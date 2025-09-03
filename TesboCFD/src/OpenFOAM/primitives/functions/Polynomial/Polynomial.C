/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial()
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(Zero),
    logActive_(false),
    logCoeff_(0)
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(std::initializer_list<scalar> coeffs)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0)
{
    if (coeffs.size() != PolySize)
    {
        FatalErrorInFunction
            << "Size mismatch: Needed " << PolySize
            << " but given " << label(coeffs.size())
            << nl << exit(FatalError);
    }

    auto iter = coeffs.begin();
    for (int i=0; i<PolySize; ++i)
    {
        this->v_[i] = *iter;
        ++iter;
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const scalar coeffs[PolySize])
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0)
{
    for (int i=0; i<PolySize; ++i)
    {
        this->v_[i] = coeffs[i];
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const UList<scalar>& coeffs)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0)
{
    if (coeffs.size() != PolySize)
    {
        FatalErrorInFunction
            << "Size mismatch: Needed " << PolySize
            << " but given " << coeffs.size()
            << nl << exit(FatalError);
    }

    for (int i = 0; i < PolySize; ++i)
    {
        this->v_[i] = coeffs[i];
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(is),
    logActive_(false),
    logCoeff_(0)
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const word& name, Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0)
{
    const word isName(is);

    if (isName != name)
    {
        FatalErrorInFunction
            << "Expected polynomial name " << name << " but read " << isName
            << nl << exit(FatalError);
    }

    is >>
        static_cast
        <
            VectorSpace<Polynomial<PolySize>, scalar, PolySize>&
        >(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int PolySize>
bool Foam::Polynomial<PolySize>::logActive() const noexcept
{
    return logActive_;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::logCoeff() const noexcept
{
    return logCoeff_;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::value(const scalar x) const
{
    scalar val = this->v_[0];

    // Avoid costly pow() in calculation
    scalar powX = x;
    for (label i=1; i<PolySize; ++i)
    {
        val += this->v_[i]*powX;
        powX *= x;
    }

    if (logActive_)
    {
        val += logCoeff_*log(x);
    }

    return val;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::derivative(const scalar x) const
{
    scalar deriv = 0;

    if (PolySize > 1)
    {
        // Avoid costly pow() in calculation
        deriv += this->v_[1];

        scalar powX = x;
        for (label i=2; i<PolySize; ++i)
        {
            deriv += i*this->v_[i]*powX;
            powX *= x;
        }
    }

    if (logActive_)
    {
        deriv += logCoeff_/x;
    }

    return deriv;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    // Avoid costly pow() in calculation
    scalar powX1 = x1;
    scalar powX2 = x2;

    scalar integ = this->v_[0]*(powX2 - powX1);
    for (label i=1; i<PolySize; ++i)
    {
        powX1 *= x1;
        powX2 *= x2;
        integ += this->v_[i]/(i + 1)*(powX2 - powX1);
    }

    if (logActive_)
    {
        integ += logCoeff_*((x2*log(x2) - x2) - (x1*log(x1) - x1));
    }

    return integ;
}


template<int PolySize>
typename Foam::Polynomial<PolySize>::intPolyType
Foam::Polynomial<PolySize>::integral(const scalar intConstant) const
{
    intPolyType newCoeffs;

    newCoeffs[0] = intConstant;
    forAll(*this, i)
    {
        newCoeffs[i+1] = this->v_[i]/(i + 1);
    }

    return newCoeffs;
}


template<int PolySize>
typename Foam::Polynomial<PolySize>::polyType
Foam::Polynomial<PolySize>::integralMinus1(const scalar intConstant) const
{
    polyType newCoeffs;

    if (this->v_[0] > VSMALL)
    {
        newCoeffs.logActive_ = true;
        newCoeffs.logCoeff_ = this->v_[0];
    }

    newCoeffs[0] = intConstant;
    for (label i=1; i<PolySize; ++i)
    {
        newCoeffs[i] = this->v_[i]/i;
    }

    return newCoeffs;
}


// ************************************************************************* //

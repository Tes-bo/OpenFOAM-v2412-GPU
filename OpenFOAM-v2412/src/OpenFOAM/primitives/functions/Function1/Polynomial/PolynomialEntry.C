/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "PolynomialEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Polynomial<Type>::checkCoefficients()
{
    if (coeffs_.empty())
    {
        FatalErrorInFunction
            << "Invalid (empty) polynomial coefficients for "
            << this->name() << nl
            << exit(FatalError);
    }

    for (const auto& coeff : coeffs_)
    {
        if (mag(coeff.second() + pTraits<Type>::one) < ROOTVSMALL)
        {
            canIntegrate_ = false;
            break;
        }
    }

    if (debug && !canIntegrate_)
    {
        WarningInFunction
            << "Polynomial " << this->name() << " cannot be integrated"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Polynomial<Type>::Polynomial
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, dict, obrPtr),
    coeffs_(),
    canIntegrate_(true)
{
    const entry* eptr = dict.findEntry(entryName, keyType::LITERAL);

    if (eptr && eptr->isStream())
    {
        // Primitive (inline) format. Eg,
        // key polynomial ((0 0) (10 1));

        ITstream& is = eptr->stream();
        if (is.peek().isWord())
        {
            is.skip();  // Discard leading 'polynomial'
        }
        is >> this->coeffs_;
        dict.checkITstream(is, entryName);
    }
    else
    {
        // Dictionary format - "coeffs" lookup. Eg,
        //
        // key { type polynomial; coeffs ((0 0) (10 1)); }

        dict.readEntry("coeffs", this->coeffs_);
    }

    this->checkCoefficients();
}


template<class Type>
Foam::Function1Types::Polynomial<Type>::Polynomial
(
    const word& entryName,
    const List<Tuple2<Type, Type>>& coeffs,
    const objectRegistry* obrPtr
)
:
    Function1<Type>(entryName, obrPtr),
    coeffs_(coeffs),
    canIntegrate_(true)
{
    this->checkCoefficients();
}


template<class Type>
Foam::Function1Types::Polynomial<Type>::Polynomial(const Polynomial& poly)
:
    Function1<Type>(poly),
    coeffs_(poly.coeffs_),
    canIntegrate_(poly.canIntegrate_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Polynomial<Type>::userTimeToTime(const Time& t)
{
    forAll(coeffs_, i)
    {
        Type value = coeffs_[i].first();
        for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
        {
            setComponent(coeffs_[i].first(), cmpt) =
                t.userTimeToTime(component(value, cmpt));
        }
    }
}


template<class Type>
Type Foam::Function1Types::Polynomial<Type>::value(const scalar x) const
{
    Type y(Zero);
    forAll(coeffs_, i)
    {
        y += cmptMultiply
        (
            coeffs_[i].first(),
            cmptPow(pTraits<Type>::one*x, coeffs_[i].second())
        );
    }

    return y;
}


template<class Type>
Type Foam::Function1Types::Polynomial<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    Type intx(Zero);

    if (canIntegrate_)
    {
        forAll(coeffs_, i)
        {
            intx += cmptMultiply
            (
                cmptDivide
                (
                    coeffs_[i].first(),
                    coeffs_[i].second() + pTraits<Type>::one
                ),
                cmptPow
                (
                    pTraits<Type>::one*x2,
                    coeffs_[i].second() + pTraits<Type>::one
                )
              - cmptPow
                (
                    pTraits<Type>::one*x1,
                    coeffs_[i].second() + pTraits<Type>::one
                )
            );
        }
    }

    return intx;
}


template<class Type>
void Foam::Function1Types::Polynomial<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << nl << indent << coeffs_;
    os.endEntry();
}


// ************************************************************************* //

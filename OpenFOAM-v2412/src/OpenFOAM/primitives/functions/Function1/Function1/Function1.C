/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "Function1.H"
#include "Time.H"
// Required by clang 5 for correct instantiation of Function1::New
#include "Constant.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::Function1
(
    const word& entryName,
    const objectRegistry* obrPtr
)
:
    function1Base(entryName, obrPtr)
{}


template<class Type>
Foam::Function1<Type>::Function1
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    function1Base(entryName, dict, obrPtr)
{}


template<class Type>
Foam::Function1<Type>::Function1(const Function1<Type>& rhs)
:
    function1Base(rhs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1<Type>::value(const scalar x) const
{
    NotImplemented;
    return Type();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::value
(
    const scalarField& x
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
Type Foam::Function1<Type>::integrate(const scalar x1, const scalar x2) const
{
    NotImplemented;
    return Type();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    NotImplemented;
    return nullptr;
}


template<class Function1Type>
Foam::tmp<Foam::Field<typename Function1Type::returnType>>
Foam::FieldFunction1<Function1Type>::value
(
    const scalarField& x
) const
{
    auto tfld = tmp<Field<Type>>::New(x.size());
    auto& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = Function1Type::value(x[i]);
    }
    return tfld;
}


template<class Function1Type>
Foam::FieldFunction1<Function1Type>::FieldFunction1
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
:
    Function1Type(entryName, dict, obrPtr)
{}


template<class Function1Type>
Foam::tmp<Foam::Field<typename Function1Type::returnType>>
Foam::FieldFunction1<Function1Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    auto tfld = tmp<Field<Type>>::New(x1.size());
    auto& fld = tfld.ref();

    forAll(x1, i)
    {
        fld[i] = Function1Type::integrate(x1[i], x2[i]);
    }

    return tfld;
}


template<class Type>
void Foam::Function1<Type>::writeEntries(Ostream& os) const
{}


template<class Type>
void Foam::Function1<Type>::writeData(Ostream& os) const
{
    os.writeKeyword(name_) << type();
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function1<Type>& rhs
)
{
    os.check(FUNCTION_NAME);

    os  << rhs.name();
    rhs.writeData(os);

    return os;
}


// ************************************************************************* //

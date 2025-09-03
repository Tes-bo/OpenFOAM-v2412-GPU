/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "PatchFunction1Expression.H"
#include "fvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::PatchExprField<Type>::PatchExprField
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    dict_(dict),  // Deep copy
    valueExpr_("expression", dict_),
    driver_(fvPatch::lookupPatch(this->patch()), dict_)
{
    if (dict.getOrDefault("debug", false))
    {
        debug |= 1;
    }

    // Basic sanity
    if (valueExpr_.empty())
    {
        FatalIOErrorInFunction(dict_)
            << "The expression was not defined!" << nl
            << exit(FatalIOError);
    }

    driver_.readDict(dict_);
}


template<class Type>
Foam::PatchFunction1Types::PatchExprField<Type>::PatchExprField
(
    const PatchExprField<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    dict_(rhs.dict_),  // Deep copy
    valueExpr_(rhs.valueExpr_),
    driver_(fvPatch::lookupPatch(this->patch()), rhs.driver_, dict_)
{}


template<class Type>
Foam::PatchFunction1Types::PatchExprField<Type>::PatchExprField
(
    const PatchExprField<Type>& rhs
)
:
    PatchExprField<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::PatchExprField<Type>::value
(
    const scalar x
) const
{
    // Expression evaluation
    driver_.clearVariables();

    driver_.setArgument(x);

    tmp<Field<Type>> tresult(driver_.evaluate<Type>(this->valueExpr_));

    DebugInfo
        << "Evaluated: " << tresult() << nl;

    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::PatchExprField<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
void Foam::PatchFunction1Types::PatchExprField<Type>::autoMap
(
    const FieldMapper& mapper
)
{
    PatchFunction1<Type>::autoMap(mapper);
}


template<class Type>
void Foam::PatchFunction1Types::PatchExprField<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{
    PatchFunction1<Type>::rmap(pf1, addr);
}


template<class Type>
void Foam::PatchFunction1Types::PatchExprField<Type>::writeData
(
    Ostream& os
) const
{
    // PatchFunction1-from-subdict so output dictionary contains
    // only the relevant entries.
    dict_.writeEntry(this->name(), os);
}


// ************************************************************************* //

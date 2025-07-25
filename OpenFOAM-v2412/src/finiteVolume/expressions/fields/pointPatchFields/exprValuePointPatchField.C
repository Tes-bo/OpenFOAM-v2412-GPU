/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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

\*---------------------------------------------------------------------------*/

#include "exprValuePointPatchField.H"
#include "pointPatchFieldMapper.H"
#include "facePointPatch.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    parent_bctype(p, iF),
    expressions::patchExprFieldBase(),
    dict_(),
    driver_
    (
        fvPatch::lookupPatch
        (
            refCast<const facePointPatch>(this->patch()).patch()
        )
    )
{}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const exprValuePointPatchField<Type>& rhs,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_
    (
        fvPatch::lookupPatch
        (
            refCast<const facePointPatch>(this->patch()).patch()
        ),
        rhs.driver_,
        dict_
    )
{}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF),  // bypass dictionary constructor
    expressions::patchExprFieldBase
    (
        dict,
        expressions::patchExprFieldBase::expectedTypes::VALUE_TYPE,
        true  // pointValue
    ),
    dict_
    (
        // Copy dictionary without "heavy" data chunks
        dictionaryContent::copyDict
        (
            dict,
            wordList(),  // allow
            wordList     // deny
            ({
                "type",  // redundant
                "value"
            })
        )
    ),
    driver_
    (
        fvPatch::lookupPatch
        (
            refCast<const facePointPatch>(this->patch()).patch()
        ),
        dict_
    )
{
    //FUTURE?
    //DeprecatedInFunction(2212)
    //    << "Use uniformFixedValue with an expression Function1 instead." << nl
    //    << "    This boundary condition will be removed in the future" << endl;

    // Require valueExpr
    if (this->valueExpr_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "The valueExpr was not defined!" << nl
            << exit(FatalIOError);
    }

    driver_.readDict(dict_);

    if (!this->readValueEntry(dict))
    {
        // Ensure field has reasonable initial values
        this->extrapolateInternal();
    }

    if (this->evalOnConstruct_)
    {
        // For potentialFoam or other solvers that don't evaluate
        this->evaluate();
    }
}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const exprValuePointPatchField<Type>& rhs,
    const DimensionedField<Type, pointMesh>& iF
)
:
    parent_bctype(rhs, iF),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_
    (
        fvPatch::lookupPatch
        (
            refCast<const facePointPatch>(this->patch()).patch()
        ),
        rhs.driver_,
        dict_
    )
{}


template<class Type>
Foam::exprValuePointPatchField<Type>::exprValuePointPatchField
(
    const exprValuePointPatchField<Type>& rhs
)
:
    parent_bctype(rhs),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_
    (
        fvPatch::lookupPatch
        (
            refCast<const facePointPatch>(this->patch()).patch()
        ),
        rhs.driver_,
        dict_
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::exprValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (debug)
    {
        InfoInFunction
            << "Value: " << this->valueExpr_ << nl
            << "Variables: ";
        driver_.writeVariableStrings(Info) << nl;
        Info<< "... updating" << endl;
    }


    // Expression evaluation
    {
        bool evalValue = (!this->valueExpr_.empty() && this->valueExpr_ != "0");


        driver_.clearVariables();

        if (evalValue)
        {
            Field<Type>::operator=
            (
                driver_.evaluate<Type>(this->valueExpr_, true)
            );
        }
        else
        {
            (*this) == Zero;
        }
    }

    this->parent_bctype::updateCoeffs();
}


template<class Type>
void Foam::exprValuePointPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    expressions::patchExprFieldBase::write(os);

    this->writeValueEntry(os);

    driver_.writeCommon(os, this->debug_ || debug);
}


// ************************************************************************* //

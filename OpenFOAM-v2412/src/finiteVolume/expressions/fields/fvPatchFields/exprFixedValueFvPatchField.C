/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2009-2018 Bernhard Gschaider
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "exprFixedValueFvPatchField.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::exprFixedValueFvPatchField<Type>::setDebug()
{
    if (expressions::patchExprFieldBase::debug_ && !debug)
    {
        debug = 1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(p, iF),
    expressions::patchExprFieldBase(),
    dict_(),
    driver_(this->patch())
{}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const exprFixedValueFvPatchField<Type>& rhs,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_(this->patch(), rhs.driver_, dict_)
{
    setDebug();
    DebugInFunction << nl;
}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue  // (ignored)
)
:
    parent_bctype(p, iF), // bypass dictionary constructor
    expressions::patchExprFieldBase
    (
        dict,
        expressions::patchExprFieldBase::expectedTypes::VALUE_TYPE
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
    driver_(this->patch(), dict_)
{
    DeprecatedInFunction(2212)
        << "Use uniformFixedValue with an expression Function1 instead." << nl
        << "    This boundary condition will be removed in the future" << endl;

    setDebug();
    DebugInFunction << nl;

    // Require valueExpr
    if (this->valueExpr_.empty())
    {
        FatalIOErrorInFunction(dict)
            << "The valueExpr was not defined!" << nl
            << exit(FatalIOError);
    }


    driver_.readDict(dict_);

    // Since we bypassed dictionary constructor
    fvPatchFieldBase::readDict(dict);

    if (!this->readValueEntry(dict))
    {
        // Ensure field has reasonable initial values
        this->extrapolateInternal();

        #ifdef FULLDEBUG
        WarningInFunction
            << "No value defined for "
            << this->internalField().name() << " on "
            << this->patch().name() << " - using patch internal field" << endl;
        #endif
    }

    if (this->evalOnConstruct_)
    {
        // For potentialFoam or other solvers that don't evaluate
        this->evaluate();
    }
}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const exprFixedValueFvPatchField<Type>& rhs
)
:
    parent_bctype(rhs),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_(this->patch(), rhs.driver_, dict_)
{
    setDebug();
    DebugInFunction << nl;
}


template<class Type>
Foam::exprFixedValueFvPatchField<Type>::exprFixedValueFvPatchField
(
    const exprFixedValueFvPatchField<Type>& rhs,
    const DimensionedField<Type, volMesh>& iF
)
:
    parent_bctype(rhs, iF),
    expressions::patchExprFieldBase(rhs),
    dict_(rhs.dict_),  // Deep copy
    driver_(this->patch(), rhs.driver_, dict_)
{
    setDebug();
    DebugInFunction << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::exprFixedValueFvPatchField<Type>::updateCoeffs()
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
            (*this) == driver_.evaluate<Type>(this->valueExpr_);
        }
        else
        {
            (*this) == Zero;
        }
    }

    this->parent_bctype::updateCoeffs();
}


template<class Type>
void Foam::exprFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    this->parent_bctype::write(os);
    expressions::patchExprFieldBase::write(os);

    driver_.writeCommon(os, this->debug_ || debug);
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
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

#include "CONSTRUCT.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::CLASS::t() const
{
    return this->db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF
)
:
    PARENT(p, iF),
    scalarData_(0),
    data_(Zero),
    fieldData_(p.size(), Zero),
    timeVsData_(),
    wordData_("wordDefault"),
    labelData_(-1),
    boolData_(false)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = Zero;
}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF,
    const dictionary& dict
)
:
    PARENT(p, iF),
    scalarData_(dict.get<scalar>("scalarData")),
    data_(dict.get<TYPE>("data")),
    fieldData_("fieldData", dict, p.size()),
    timeVsData_(Function1<TYPE>::New("timeVsData", dict, &this->db())),
    wordData_(dict.getOrDefault<word>("wordName", "wordDefault")),
    labelData_(-1),
    boolData_(false)
{
    this->refGrad() = Zero;
    this->valueFraction() = Zero;

    this->refValue() = FIELD("fieldData", dict, p.size());
    FVPATCHF::operator=(this->refValue());

    PARENT::evaluate();

    /*
    //Initialise with the value entry if evaluation is not possible
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
    this->refValue() = *this;
    */
}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf,
    const fvPatch& p,
    const DimensionedField<TYPE, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    PARENT(ptf, p, iF, mapper),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_, mapper),
    timeVsData_(ptf.timeVsData_.clone()),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf
)
:
    PARENT(ptf),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_.clone()),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


template<class Type>
Foam::CLASS::
CONSTRUCT
(
    const CLASS& ptf,
    const DimensionedField<TYPE, volMesh>& iF
)
:
    PARENT(ptf, iF),
    scalarData_(ptf.scalarData_),
    data_(ptf.data_),
    fieldData_(ptf.fieldData_),
    timeVsData_(ptf.timeVsData_.clone()),
    wordData_(ptf.wordData_),
    labelData_(-1),
    boolData_(ptf.boolData_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::CLASS::autoMap
(
    const fvPatchFieldMapper& m
)
{
    PARENT::autoMap(m);
    fieldData_.autoMap(m);
}


template<class Type>
void Foam::CLASS::rmap
(
    const FVPATCHF& ptf,
    const labelList& addr
)
{
    PARENT::rmap(ptf, addr);

    const CLASS& tiptf =
        refCast<const CLASS>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}


template<class Type>
void Foam::CLASS::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    PARENT::operator==
    (
        data_
      + fieldData_
      + scalarData_*timeVsData_->value(t())
    );

    const scalarField& phip =
        this->patch().template lookupPatchField<surfaceScalarField>("phi");

    this->valueFraction() = neg(phip);

    PARENT::updateCoeffs();
}


template<class Type>
void Foam::CLASS::write
(
    Ostream& os
) const
{
    FVPATCHF::write(os);
    os.writeEntry("scalarData", scalarData_);
    os.writeEntry("data", data_);
    fieldData_.writeEntry("fieldData", os);
    timeVsData_->writeData(os);
    os.writeEntry("wordData", wordData_);
    FVPATCHF::writeValueEntry(os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        FVPATCHF,
        CLASS
    );
}

// ************************************************************************* //

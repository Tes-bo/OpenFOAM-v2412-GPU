/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "edgeNormalFixedValueFaPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "areaFields.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(p, iF),
    refValue_(p.size(), Zero)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const scalar refValue
)
:
    fixedValueFaPatchVectorField(p, iF),
    refValue_(p.size(), refValue)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& ptf,
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchVectorField(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const faPatch& p,
    const DimensionedField<vector, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchVectorField(p, iF, dict, IOobjectOption::NO_READ),
    refValue_("refValue", dict, p.size())
{
    tmp<vectorField> tvalues(refValue_*patch().edgeNormals());

    faPatchVectorField::operator=(tvalues);
}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& pivpvf
)
:
    fixedValueFaPatchVectorField(pivpvf),
    refValue_(pivpvf.refValue_)
{}


Foam::edgeNormalFixedValueFaPatchVectorField::
edgeNormalFixedValueFaPatchVectorField
(
    const edgeNormalFixedValueFaPatchVectorField& pivpvf,
    const DimensionedField<vector, areaMesh>& iF
)
:
    fixedValueFaPatchVectorField(pivpvf, iF),
    refValue_(pivpvf.refValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::edgeNormalFixedValueFaPatchVectorField::autoMap
(
    const faPatchFieldMapper& m
)
{
    fixedValueFaPatchVectorField::autoMap(m);
    refValue_.autoMap(m);
}


void Foam::edgeNormalFixedValueFaPatchVectorField::rmap
(
    const faPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFaPatchVectorField::rmap(ptf, addr);

    const edgeNormalFixedValueFaPatchVectorField& tiptf =
        refCast<const edgeNormalFixedValueFaPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::edgeNormalFixedValueFaPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    tmp<vectorField> tvalues(refValue_*patch().edgeNormals());

    faPatchVectorField::operator=(tvalues);
    faPatchVectorField::updateCoeffs();
}


void Foam::edgeNormalFixedValueFaPatchVectorField::write(Ostream& os) const
{
    fixedValueFaPatchField<vector>::write(os);
    refValue_.writeEntry("refValue", os);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeFaPatchTypeField
    (
        faPatchVectorField,
        edgeNormalFixedValueFaPatchVectorField
    );
}


// ************************************************************************* //

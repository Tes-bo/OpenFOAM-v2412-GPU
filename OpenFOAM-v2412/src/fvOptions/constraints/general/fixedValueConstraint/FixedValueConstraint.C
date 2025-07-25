/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "FixedValueConstraint.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::FixedValueConstraint<Type>::FixedValueConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::fv::FixedValueConstraint<Type>::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        const dictionary& fieldValuesDict = coeffs_.subDict("fieldValues");

        label count = fieldValuesDict.size();

        fieldNames_.resize(count);
        fieldValues_.resize(count);

        fv::option::resetApplied();

        count = 0;
        for (const entry& dEntry : fieldValuesDict)
        {
            fieldNames_[count] = dEntry.keyword();
            dEntry.readEntry(fieldValues_[count]);

            ++count;
        }

        return true;
    }

    return false;
}


template<class Type>
void Foam::fv::FixedValueConstraint<Type>::constrain
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    DebugInfo
        << "FixedValueConstraint<"
        << pTraits<Type>::typeName
        << ">::constrain for source " << name_ << endl;

    eqn.setValues(cells_, fieldValues_[fieldi]);
}


// ************************************************************************* //

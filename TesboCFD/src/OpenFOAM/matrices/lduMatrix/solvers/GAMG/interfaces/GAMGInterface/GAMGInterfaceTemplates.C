/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "GAMGInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::GAMGInterface::interfaceInternalField
(
    const UList<Type>& iF
) const
{
    auto tresult = tmp<Field<Type>>::New(size());
    interfaceInternalField(iF, tresult.ref());
    return tresult;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::GAMGInterface::interfaceInternalField
(
    const UList<Type>& iF,
    const labelUList& faceCells
) const
{
    auto tresult = tmp<Field<Type>>::New(faceCells.size());
    auto& result = tresult.ref();

    forAll(result, elemi)
    {
        result[elemi] = iF[faceCells[elemi]];
    }
    return tresult;
}


template<class Type>
void Foam::GAMGInterface::interfaceInternalField
(
    const UList<Type>& iF,
    List<Type>& result
) const
{
    result.resize(size());

    forAll(result, elemi)
    {
        result[elemi] = iF[faceCells_[elemi]];
    }
}


// ************************************************************************* //

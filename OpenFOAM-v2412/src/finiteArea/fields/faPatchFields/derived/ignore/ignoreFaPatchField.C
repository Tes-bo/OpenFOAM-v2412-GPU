/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>

\*---------------------------------------------------------------------------*/

#include "ignoreFaPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::ignoreFaPatchField<Type>::ignoreFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    zeroGradientFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::ignoreFaPatchField<Type>::ignoreFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFaPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::ignoreFaPatchField<Type>::ignoreFaPatchField
(
    const ignoreFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    zeroGradientFaPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::ignoreFaPatchField<Type>::ignoreFaPatchField
(
    const ignoreFaPatchField<Type>& ptf
)
:
    zeroGradientFaPatchField<Type>(ptf)
{}


template<class Type>
Foam::ignoreFaPatchField<Type>::ignoreFaPatchField
(
    const ignoreFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    zeroGradientFaPatchField<Type>(ptf, iF)
{}


// ************************************************************************* //

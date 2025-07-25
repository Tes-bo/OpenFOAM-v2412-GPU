/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "OStringStream.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::string Foam::exprTools::zeroValue()
{
    OStringStream buf;

    buf << pTraits<Type>::typeName << '(';
    for (direction cmpt=0;  cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        if (cmpt) buf << ',';
        buf << 0;
    }
    buf << ')';

    return buf.str();
}


template<class Type>
Foam::string Foam::exprTools::toString
(
    const Type& data,
    const word& prefix
)
{
    OStringStream buf;

    buf << prefix << '(';
    for (direction cmpt=0;  cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        if (cmpt) buf << ',';
        buf << component(data, cmpt);
    }
    buf << ')';

    return buf.str();
}


template<class Type>
Foam::string Foam::exprTools::toString(const Type& data)
{
    return toString<Type>(data, pTraits<Type>::typeName);
}


template<class Type>
Foam::string Foam::exprTools::toString(ITstream& is)
{
    Type data(Zero);
    is >> data;

    return toString<Type>(data, pTraits<Type>::typeName);
}


// ************************************************************************* //

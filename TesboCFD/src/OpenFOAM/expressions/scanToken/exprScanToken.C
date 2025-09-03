/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "exprScanToken.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::expressions::scanToken Foam::expressions::scanToken::null()
{
    scanToken tok;
    tok.type_ = tokenType::LABEL;
    tok.labelValue = 0;

    return tok;
}


void Foam::expressions::scanToken::destroy()
{
    if (type_ == tokenType::VECTOR)
    {
        delete vectorPtr;
        vectorPtr = nullptr;
    }
    else if (type_ == tokenType::WORD)
    {
        delete wordPtr;
        wordPtr = nullptr;
    }
}


void Foam::expressions::scanToken::setLabel(label val)
{
    type_ = tokenType::LABEL;
    labelValue = val;
}


void Foam::expressions::scanToken::setScalar(scalar val)
{
    type_ = tokenType::SCALAR;
    scalarValue = val;
}


void Foam::expressions::scanToken::setVector(scalar x, scalar y, scalar z)
{
    type_ = tokenType::VECTOR;
    vectorPtr = new Foam::vector(x, y, z);
}


void Foam::expressions::scanToken::setVector(const vector& val)
{
    type_ = tokenType::VECTOR;
    vectorPtr = new Foam::vector(val);
}


void Foam::expressions::scanToken::setWord(const word& val)
{
    type_ = tokenType::WORD;
    wordPtr = new Foam::word(val);
}


// ************************************************************************* //

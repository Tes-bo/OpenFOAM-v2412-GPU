/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "stringOpsEvaluate.H"
#include "stringOps.H"
#include "StringStream.H"
#include "fieldExprDriver.H"
#include "error.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::string Foam::stringOps::evaluate
(
    label fieldWidth,
    const std::string& str,
    size_t pos,
    size_t len
)
{
    if (fieldWidth < 1)
    {
        fieldWidth = 1;
    }

    /// InfoErr<< "Evaluate " << str.substr(pos, len) << nl;

    const auto trimPoints = stringOps::findTrim(str, pos, len);

    pos = trimPoints.first;
    len = (trimPoints.second - trimPoints.first);

    if (!len)
    {
        return string();
    }

    /// InfoErr<< "Evaluate " << str.substr(pos, len) << nl;

    expressions::exprResult result;
    {
        expressions::fieldExprDriver driver(fieldWidth);
        driver.parse(str, pos, len);
        result = std::move(driver.result());
    }

    if (!result.hasValue() || !result.size())
    {
        InfoErr
            << "Failed evaluation: "
            << str.substr(pos, len) << nl;

        return string();
    }

    OStringStream os;
    if (result.size() <= 1)
    {
        result.writeValue(os);
    }
    else
    {
        result.writeField(os);
    }

    return os.str();
}


Foam::string Foam::stringOps::evaluate
(
    const std::string& str,
    size_t pos,
    size_t len
)
{
    return stringOps::evaluate(1, str, pos, len);
}


// ************************************************************************* //

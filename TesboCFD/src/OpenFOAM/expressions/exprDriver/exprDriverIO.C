/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "exprDriver.H"
#include "exprTools.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//! \cond file_scope
//  Write list as single or multiple entries - see exprTools::getList()
static void writeList
(
    Ostream& os,
    const UList<expressions::exprString>& list
)
{
    if (list.size() == 1)
    {
        os << list[0];
    }
    else
    {
        os << token::BEGIN_LIST;

        if (!list.empty())
        {
            os << nl;

            for (const expressions::exprString& str : list)
            {
                os << str << nl;
            }
        }
        os << token::END_LIST;
    }
}

//! \endcond

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::expressions::exprString
Foam::expressions::exprDriver::readExpression
(
    const word& name,
    const dictionary& dict
)
{
    return expressions::exprString(name, dict);
}


Foam::expressions::exprString
Foam::expressions::exprDriver::readExpression
(
    const word& name
)
{
    return expressions::exprString(name, dict());
}


Foam::List<Foam::expressions::exprString>
Foam::expressions::exprDriver::readVariableStrings
(
    const dictionary& dict,
    const word& keyword,
    bool mandatory
)
{
    return exprTools::getList(keyword, dict, mandatory);
}


Foam::label Foam::expressions::exprDriver::setVariableStrings
(
    const dictionary& dict,
    bool mandatory
)
{
    variableStrings_ = readVariableStrings(dict, "variables", mandatory);

    return variableStrings_.size();
}


Foam::Ostream& Foam::expressions::exprDriver::writeVariableStrings
(
    Ostream& os,
    const word& keyword
) const
{
    if (keyword.size())
    {
        os.writeKeyword(keyword);
    }

    writeList(os, variableStrings_);

    if (keyword.size())
    {
        os.endEntry();
    }

    return os;
}


// ************************************************************************* //

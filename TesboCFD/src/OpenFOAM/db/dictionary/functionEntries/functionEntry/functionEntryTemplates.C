/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include "functionEntry.H"
#include "ISstream.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class StringType>
Foam::List<StringType>
Foam::functionEntry::readStringList(Istream& is)
{
    List<StringType> list;

    ISstream& iss = dynamic_cast<ISstream&>(is);
    token firstToken(iss);

    const bool isStringTok = firstToken.isStringType();

    iss.putBack(firstToken);

    if (isStringTok)
    {
        // The first token appears viable as non-list
        // - treated like list with a single entry

        list.resize(1);
        iss >> list.front();
    }
    else
    {
        iss >> list;
    }

    return list;
}


// ************************************************************************* //

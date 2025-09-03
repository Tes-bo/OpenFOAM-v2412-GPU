/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "UILList.H"
#include "Ostream.H"
#include "token.H"

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class LListBase, class T>
Foam::Ostream& Foam::UILList<LListBase, T>::writeList
(
    Ostream& os,
    const label shortLen
) const
{
    // NB: no binary, contiguous output

    const label len = this->size();

    if
    (
        (len <= 1 || !shortLen)
     || (len <= shortLen)
    )
    {
        // Size and start delimiter
        os << len << token::BEGIN_LIST;

        // Contents
        bool space = false;
        for (const T& val : *this)
        {
            if (space) os << token::SPACE;
            space = true;
            os << val;
        }

        // End delimiter
        os << token::END_LIST;
    }
    else
    {
        // Size and start delimiter
        os << nl << len << nl << token::BEGIN_LIST << nl;

        // Contents
        for (const T& val : *this)
        {
            os << val << nl;
        }

        // End delimiter
        os << token::END_LIST;
    }

    os.check(FUNCTION_NAME);
    return os;

}


template<class LListBase, class T>
Foam::Ostream& Foam::operator<<(Ostream& os, const UILList<LListBase, T>& lst)
{
    return lst.writeList(os, -1);
}


// ************************************************************************* //

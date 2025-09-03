/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::blockMeshTools::read
(
    Istream& is,
    List<T>& list,
    const dictionary& dict
)
{
    token tok(is);

    if (tok.isLabel())
    {
        // Label: should be int(..)

        const label len = tok.labelToken();

        // Resize to length required
        list.resize_nocopy(len);

        // Begin of contents marker
        const char delimiter = is.readBeginList("List");

        if (len)
        {
            if (delimiter == token::BEGIN_LIST)
            {
                for (label i=0; i<len; ++i)
                {
                    read(is, list[i], dict);
                }
            }
        }

        // End of contents marker
        is.readEndList("List");
    }
    else if (tok.isPunctuation(token::BEGIN_LIST))
    {
        // "(...)" : read with DynamicList and transfer contents

        DynamicList<T> elements(std::move(list));
        elements.clear();  // Reset addressing only

        is >> tok;
        is.fatalCheck(FUNCTION_NAME);

        while (!tok.isPunctuation(token::END_LIST))
        {
            is.putBack(tok);

            read(is, elements.emplace_back(), dict);

            is >> tok;
            is.fatalCheck(FUNCTION_NAME);
        }

        // Transfer back to regular list
        list = std::move(elements);
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << tok.info() << nl
            << exit(FatalIOError);
    }
}


template<class T>
Foam::List<T> Foam::blockMeshTools::read
(
    Istream& is,
    const dictionary& dict
)
{
    List<T> list;
    read(is, list, dict);
    return list;
}


// ************************************************************************* //

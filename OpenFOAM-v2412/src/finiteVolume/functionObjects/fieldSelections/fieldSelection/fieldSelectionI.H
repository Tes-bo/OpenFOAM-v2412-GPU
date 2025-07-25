/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

inline Foam::HashSet<Foam::wordRe>
Foam::functionObjects::fieldSelection::filters() const
{
    HashSet<wordRe> values(2*this->size());
    for (const fieldInfo& fi : *this)
    {
        values.insert(fi.name());
    }

    return values;
}


inline Foam::wordList
Foam::functionObjects::fieldSelection::selectionNames() const
{
    DynamicList<word> values(selection_.size());

    for (const fieldInfo& fi : selection_)
    {
        values.push_uniq(fi.name());
    }

    wordList names(std::move(values));
    Foam::sort(names);  // Globally consistent order

    return names;
}


// ************************************************************************* //

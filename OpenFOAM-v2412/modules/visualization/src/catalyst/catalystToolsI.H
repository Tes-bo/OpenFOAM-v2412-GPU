/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

inline vtkCPDataDescription*
Foam::catalyst::dataQuery::description() const
{
    return descrip_;
}


inline void Foam::catalyst::dataQuery::clear()
{
    names_.clear();
    fields_.clear();
    active_.clear();
}


inline bool Foam::catalyst::dataQuery::empty() const noexcept
{
    return active_.empty();
}


inline Foam::label Foam::catalyst::dataQuery::size() const noexcept
{
    return active_.size();
}


inline bool Foam::catalyst::dataQuery::contains(const word& channel) const
{
    return fields_.contains(channel);
}


inline Foam::List<Foam::label> Foam::catalyst::dataQuery::active() const
{
    // The active indices
    return active_.sortedToc();
}


inline Foam::wordList Foam::catalyst::dataQuery::names() const
{
    // Subset of known names with the active indices
    return wordList(names_, active_.sortedToc());
}


inline const Foam::wordHashSet&
Foam::catalyst::dataQuery::fields(const word& channel) const
{
    const auto iter = fields_.cfind(channel);
    if (iter.good())
    {
        return *iter;
    }

    return emptyWordHash_;
}


inline void Foam::catalyst::dataQuery::clear(const word& channel)
{
    fields_.erase(channel);
    active_.erase(names_.find(channel));
}


// ************************************************************************* //

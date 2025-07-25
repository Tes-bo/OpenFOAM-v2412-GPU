/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::DirLister::DirLister
(
    const fileName& directory,
    const fileName::Type requestedType,
    const bool filtergz,
    const bool followLink
)
:
    dirName_(directory),
    requestedType_(requestedType),
    followLink_(followLink),
    stripgz_(filtergz)
{}


inline Foam::DirLister Foam::DirLister::dirs
(
    const fileName& directory,
    const bool followLink
)
{
    return DirLister(directory, fileName::DIRECTORY, false, followLink);
}


inline Foam::DirLister Foam::DirLister::files
(
    const fileName& directory,
    const bool filtergz,
    const bool followLink
)
{
    return DirLister(directory, fileName::FILE, filtergz, followLink);
}


template<class UnaryPredicate>
inline Foam::DirLister::Filtered<UnaryPredicate>::Filtered
(
    const Foam::DirLister& lister,
    const UnaryPredicate& pred,
    const bool prune
)
:
    DirLister(lister),
    pred_(pred),
    prune_(prune)
{}


template<class UnaryPredicate>
inline Foam::DirLister::Filtered<UnaryPredicate>
Foam::DirLister::where(const UnaryPredicate& pred, const bool prune) const
{
    return DirLister::Filtered<UnaryPredicate>(*this, pred, prune);
}


// * * * * * * * * * * * * * * * Const Iterator  * * * * * * * * * * * * * * //

inline Foam::DirLister::const_iterator::const_iterator()
:
    lister_(nullptr),
    dirptr_(nullptr),
    name_()
{}


inline Foam::DirLister::const_iterator::const_iterator
(
    const Foam::DirLister* lister
)
:
    lister_(lister),
    dirptr_(nullptr),
    name_()
{
    if (lister_ && open(lister_->dirName()))
    {
        next();  // increment to first entry
    }
}


inline Foam::DirLister::const_iterator::~const_iterator()
{
    close();
}


Foam::DirLister::const_iterator&
Foam::DirLister::const_iterator::operator++()
{
    next();
    return *this;
}


inline bool Foam::DirLister::const_iterator::operator==
(
    const const_iterator& iter
)
const
{
    // Directory entries are unique, so just compare the names
    return (name_ == iter.name_);
}


inline bool Foam::DirLister::const_iterator::operator!=
(
    const const_iterator& iter
)
const
{
    return name_ != iter.name_;
}


inline Foam::DirLister::const_iterator Foam::DirLister::cbegin() const
{
    return const_iterator(this);
}


inline Foam::DirLister::const_iterator Foam::DirLister::cend() const
{
    return const_iterator();
}


inline Foam::DirLister::const_iterator Foam::DirLister::begin() const
{
    return cbegin();
}


inline Foam::DirLister::const_iterator Foam::DirLister::end() const
{
    return cend();
}


// ************************************************************************* //

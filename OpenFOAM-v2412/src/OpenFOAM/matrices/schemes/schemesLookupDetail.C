/*---------------------------------------------------------------------------*\
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

#include "schemesLookup.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::schemesLookup::lookupDetail::lookupDetail
(
    const word& dictName,
    const fileName& parentDictPath
)
:
    name_(dictName),
    dict_(),
    default_()
{
    if (parentDictPath.empty())
    {
        dict_.name() = name_;
    }
    else if (name_.empty())
    {
        dict_.name() = parentDictPath;
        name_ = dict_.dictName();
    }
    else
    {
        dict_.name() = fileName::concat(parentDictPath, name_, '/');
    }
    default_.name() = fileName::concat(dict_.name(), "default", '/');
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::schemesLookup::lookupDetail::clear()
{
    dict_.clear();
    default_.clear();
}


Foam::ITstream& Foam::schemesLookup::lookupDetail::fallback() const
{
    ITstream& is = const_cast<ITstream&>(default_);
    is.rewind();
    return is;
}


Foam::ITstream& Foam::schemesLookup::lookupDetail::lookup
(
    const word& name
) const
{
    if (dict_.found(name) || default_.empty())
    {
        // Fails when 'name' is not found and no default is available,
        // which provides some traceback error messages
        return dict_.lookup(name);
    }

    return fallback();
}


void Foam::schemesLookup::lookupDetail::populate
(
    const dictionary& dict,
    const word& defaultName,
    const bool mandatory
)
{
    if (mandatory || dict.found(name_))
    {
        // Fails when 'name' is not found but it is mandatory,
        // which provides some traceback error messages
        dict_ = dict.subDict(name_);
    }
    else if (!defaultName.empty() && !dict_.found("default"))
    {
        dict_.add("default", defaultName);
    }

    // Clear or set the default stream
    if
    (
        !dict_.found("default")
     || dict_.lookup("default").peek() == "none"
    )
    {
        default_.clear();
        default_.rewind();
    }
    else
    {
        default_ = dict_.lookup("default");
    }
}


void Foam::schemesLookup::lookupDetail::writeEntry(Ostream& os) const
{
    dict_.writeEntry(os);
}


void Foam::schemesLookup::lookupDetail::writeEntryOptional(Ostream& os) const
{
    if (!dict_.empty())
    {
        dict_.writeEntry(os);
    }
}


// ************************************************************************* //

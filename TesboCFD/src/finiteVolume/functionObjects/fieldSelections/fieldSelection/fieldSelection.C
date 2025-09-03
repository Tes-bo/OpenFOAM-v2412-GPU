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

#include "fieldSelection.H"
#include "objectRegistry.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldSelection::fieldSelection
(
    const objectRegistry& obr,
    const bool includeComponents
)
:
    obr_(obr),
    includeComponents_(includeComponents)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container>
bool Foam::functionObjects::fieldSelection::resetFieldFiltersImpl
(
    const Container& names
)
{
    static std::string cmptStr(".component(");
    static std::string::size_type len(cmptStr.size());

    DynamicList<fieldInfo> nameAndComponent(names.size());

    for (const wordRe& name : names)
    {
        const auto n = name.find(cmptStr);
        if (n != std::string::npos)
        {
            // Field should be written <field>.component(i)

            if (!includeComponents_)
            {
                FatalErrorInFunction
                    << "Component specification not allowed for " << name
                    << exit(FatalError);
            }

            if (name.isPattern())
            {
                FatalErrorInFunction
                    << "Cannot use \".component option\" in combination with "
                    << "wildcards for " << name
                    << exit(FatalError);
            }

            const word baseName(name.substr(0, n));

            // Extract the component - number between ()'s
            const auto closei = name.find(')', n);

            if (closei == std::string::npos)
            {
                FatalErrorInFunction
                    << "Invalid field component specification for "
                    << name << nl
                    << ". Field should be expressed as <field>.component(i)"
                    << exit(FatalError);
            }

            const auto cmptWidth = (closei - n - len);

            label component
            (
                readLabel(name.substr(n+len, cmptWidth))
            );

            nameAndComponent.emplace_back(wordRe(baseName), component);
        }
        else
        {
            nameAndComponent.emplace_back(name);
        }
    }

    this->transfer(nameAndComponent);

    return true;
}


bool Foam::functionObjects::fieldSelection::resetFieldFilters
(
    const HashSet<wordRe>& names
)
{
    return resetFieldFiltersImpl(names);
}


bool Foam::functionObjects::fieldSelection::resetFieldFilters
(
    const wordRe& name
)
{
    List<wordRe> names(1, name);

    return resetFieldFiltersImpl(names);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldSelection::read(const dictionary& dict)
{
    HashSet<wordRe> fields;
    dict.readEntry("fields", fields);

    return resetFieldFilters(fields);
}


bool Foam::functionObjects::fieldSelection::containsPattern() const
{
    for (const fieldInfo& fi : *this)
    {
        if (fi.name().isPattern())
        {
            return true;
        }
    }

    return false;
}


void Foam::functionObjects::fieldSelection::clearSelection()
{
    selection_.clear();
}


bool Foam::functionObjects::fieldSelection::updateSelection()
{
    return false;
}


bool Foam::functionObjects::fieldSelection::checkSelection()
{
    bool ok = true;
    for (const fieldInfo& fi : *this)
    {
        if (!fi.found())
        {
            WarningInFunction
                << "Field " << fi.name() << " not found"
                << endl;

            ok = false;
        }
    }

    return ok;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "mapDistributeBase.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// The maps (labelListList) are not human-modifiable but if we need to
// inspect them in ASCII, it is much more convenient if each sub-list
// is flattened on a single line.
static Ostream& printMaps(Ostream& os, const labelListList& maps)
{
    if (os.format() == IOstreamOption::BINARY || maps.empty())
    {
        os  << maps;
    }
    else
    {
        os  << nl << maps.size() << nl
            << token::BEGIN_LIST << nl;

        // Compact single-line output for each labelList
        for (const labelList& map : maps)
        {
            map.writeList(os) << nl;
        }
        os  << token::END_LIST;
    }

    return os;
}


static void writeMaps(Ostream& os, const word& key, const labelListList& maps)
{
    if (os.format() == IOstreamOption::BINARY || maps.empty())
    {
        os.writeEntry(key, maps);
    }
    else
    {
        os  << indent << key;
        printMaps(os, maps);
        os.endEntry();
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistributeBase::mapDistributeBase
(
    const dictionary& dict,
    const label comm
)
:
    mapDistributeBase(comm)
{
    mapDistributeBase::readDict(dict);
}


Foam::mapDistributeBase::mapDistributeBase(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::mapDistributeBase::readDict(const dictionary& dict)
{
    constructSize_ = dict.get<label>("constructSize");

    // The subMap
    {
        const dictionary& subdict = dict.subDict("subMap");

        subdict.readEntry("flip", subHasFlip_);
        subdict.readEntry("maps", subMap_);
    }

    // The constructMap
    {
        const dictionary& subdict = dict.subDict("constructMap");

        subdict.readEntry("flip", constructHasFlip_);
        subdict.readEntry("maps", constructMap_);
    }
}


void Foam::mapDistributeBase::writeEntries(Ostream& os) const
{
    os.writeEntry("constructSize", constructSize_);

    os << nl;
    os.beginBlock("subMap");
    os.writeEntry("flip", subHasFlip_);
    writeMaps(os, "maps", subMap_);
    os.endBlock();

    os << nl;
    os.beginBlock("constructMap");
    os.writeEntry("flip", constructHasFlip_);
    writeMaps(os, "maps", constructMap_);
    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, mapDistributeBase& map)
{
    is.fatalCheck(FUNCTION_NAME);

    is  >> map.constructSize_
        >> map.subMap_ >> map.constructMap_
        >> map.subHasFlip_ >> map.constructHasFlip_
        >> map.comm_;

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const mapDistributeBase& map)
{
    os  << map.constructSize_ << token::NL;

    printMaps(os, map.subMap_) << token::NL;
    printMaps(os, map.constructMap_) << token::NL;

    os  << map.subHasFlip_ << token::SPACE
        << map.constructHasFlip_ << token::SPACE
        << map.comm_ << token::NL;

    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<mapDistributeBase>& iproxy
)
{
    const auto& map = *iproxy;

    // Output as compact pseudo dictionary entries

    os.writeEntry("constructSize", map.constructSize());

    os  << indent << "local  { flip " << map.subHasFlip()
        << "; sizes ";
    map.subMapSizes().writeList(os) << "; }" << nl;

    os  << indent << "remote { flip " << map.constructHasFlip()
        << "; sizes ";
    map.constructMapSizes().writeList(os) << "; }" << nl;

    return os;
}


// ************************************************************************* //

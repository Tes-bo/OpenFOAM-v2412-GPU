/*---------------------------------------------------------------------------*\
  =========                |
  \      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \    /   O peration     |
    \  /    A nd           | www.openfoam.com
     \/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "substitutionModel.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(substitutionModel, 0);
    defineRunTimeSelectionTable(substitutionModel, dictionary);
}

const Foam::word Foam::substitutionModel::KEY_BEGIN = "{{";
const Foam::word Foam::substitutionModel::KEY_END = "}}";
Foam::HashTable<Foam::string> Foam::substitutionModel::builtin_;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::substitutionModel::keyify(const word& w)
{
    return KEY_BEGIN + w + KEY_END;
}


Foam::word Foam::substitutionModel::cleanKey(const string& str)
{
    return stringOps::upper(stringOps::trim(str));
};


Foam::wordList Foam::substitutionModel::getKeys(string& buffer)
{
    const label lBegin = KEY_BEGIN.length();
    const label lEnd = KEY_END.length();

    wordHashSet keys;

    size_t pos0 = 0;
    size_t pos = 0;
    string cleanedBuffer = "";
    while (((pos = buffer.find(KEY_BEGIN, pos)) != string::npos))
    {
        cleanedBuffer += buffer.substr(pos0, pos-pos0);

        size_t posEnd = buffer.find(KEY_END, pos);

        if (posEnd != string::npos)
        {
            const word k(cleanKey(buffer.substr(pos+lBegin, posEnd-pos-lEnd)));
            keys.insert(k);
            cleanedBuffer += keyify(k);
        }

        pos = posEnd + lEnd;
        pos0 = pos;
    }

    cleanedBuffer += buffer.substr(pos0, buffer.length() - pos0);
    buffer = cleanedBuffer;

    return keys.toc();
}


void Foam::substitutionModel::addBuiltinStr
(
    const word& key,
    const string& value
)
{
    builtin_.insert(cleanKey(key), value.c_str());
}


bool Foam::substitutionModel::containsBuiltin(const word& key)
{
    return builtin_.contains(key);
}


void Foam::substitutionModel::setBuiltinStr
(
    const word& key,
    const string& value
)
{
    builtin_.set(cleanKey(key), value.c_str());
}


bool Foam::substitutionModel::replaceBuiltin(const word& key, string& str)
{
    if (builtin_.found(key))
    {
        str.replaceAll(keyify(key), builtin_[key].c_str());
        return true;
    }

    return false;
}


bool Foam::substitutionModel::replaceBuiltin(string& str)
{
    const string str0 = str;

    // Quick exit if there are no keys in the string
    if (str.find(KEY_BEGIN) == string::npos) return false;

    forAllConstIters(builtin_, iter)
    {
        str.replaceAll(keyify(iter.key()), iter.val().c_str());
    }

    return str != str0;
}


void Foam::substitutionModel::writeBuiltins(Ostream& os)
{
    for (const auto& iter : builtin_.csorted())
    {
        os  << keyify(iter.key()).c_str() << " : " << iter.val() << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::substitutionModel::substitutionModel
(
    const dictionary& dict,
    const Time& time
)
:
    dict_(dict),
    time_(time)
{}


// ************************************************************************* //
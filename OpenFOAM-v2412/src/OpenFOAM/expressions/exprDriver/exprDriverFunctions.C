/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "exprDriver.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Check for acceptable Function1 keywords in the given dictionary,
// with special handling to avoid accidental inclusion of coeffs
// dictionaries etc
static wordHashSet getAcceptableFunctionKeys
(
    const dictionary* dictPtr,
    const bool acceptPrimitiveEntry,
    const bool report = false
)
{
    wordHashSet acceptKeys;
    wordHashSet rejectKeys;

    if (!dictPtr)
    {
        return acceptKeys;
    }

    const dictionary& dict = *dictPtr;

    acceptKeys.reserve(dict.size());
    rejectKeys.reserve(dict.size());

    for (const entry& dEntry : dict)
    {
        const keyType& kw = dEntry.keyword();

        bool ok = true;

        if (kw.isPattern())
        {
            ok = false;
        }
        else if (dEntry.isDict())
        {
            // Dictionary entry - require "type", which should eliminate
            // any *Coeffs dictionaries

            ok = dEntry.dict().found("type", keyType::LITERAL);
        }
        else
        {
            // Primitive entry. Trust that it is okay?
            ok = acceptPrimitiveEntry;
        }

        if (ok)
        {
            acceptKeys.insert(kw);
        }
        else
        {
            rejectKeys.insert(kw);
        }
    }

    if (report && rejectKeys.size())
    {
        InfoInFunction
            << "Dropped invalid/redundant entries: "
            << flatOutput(rejectKeys.sortedToc()) << nl;
    }

    return acceptKeys;
}


// Read and reset Function1 for given dictionary.
// Uses getAcceptableFunctionKeys
template<class Type>
static void resetFuncsImpl
(
    const word& subDictName,
    const dictionary& topDict,
    HashTable<refPtr<Function1<Type>>>& tbl,
    const objectRegistry* obrPtr
)
{
    tbl.clear();

    const dictionary* dictPtr =
        topDict.findDict(subDictName, keyType::LITERAL);

    if (!dictPtr)
    {
        return;
    }

    wordHashSet acceptKeys
    (
        getAcceptableFunctionKeys
        (
            dictPtr,
            true  // Accept primitive entries, hope for the best
        )
    );

    const dictionary& dict = *dictPtr;

    for (const word& entryName : acceptKeys)
    {
        // From autoPtr -> refPtr
        refPtr<Function1<Type>> func
        (
            Function1<Type>::New(entryName, dict, obrPtr)
        );

        if (func)
        {
            tbl.insert(entryName, std::move(func));
        }
    }
}


// Write out entries, if they originated from dictionary
template<class Type>
static void writeFuncsImpl
(
    Ostream& os,
    const word& subDictName,
    const dictionary& topDict,
    const HashTable<refPtr<Function1<Type>>>& tbl
)
{
    const dictionary* dictPtr =
        topDict.findDict(subDictName, keyType::LITERAL);

    if (!dictPtr || tbl.empty())
    {
        return;
    }

    label nwrote = 0;

    const dictionary& dict = *dictPtr;

    for (const entry& dEntry : dict)
    {
        const word& entryName = dEntry.keyword();

        const auto iter = tbl.cfind(entryName);

        if (!iter.good())
        {
            continue;
        }

        const auto& funcPtr = iter.val();

        if (funcPtr)
        {
            if (!nwrote++)
            {
                os.beginBlock(subDictName);
            }

            (*funcPtr).writeData(os);
        }
    }

    if (nwrote)
    {
        os.endBlock();
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::exprDriver::resetFunctions
(
    const dictionary& dict
)
{
    resetFuncsImpl<scalar>("functions<scalar>", dict_, scalarFuncs_, obrPtr_);
    resetFuncsImpl<vector>("functions<vector>", dict_, vectorFuncs_, obrPtr_);

    if (debug)
    {
        writeFunctions(InfoInFunction);
    }
}


void Foam::expressions::exprDriver::writeFunctions(Ostream& os) const
{
    writeFuncsImpl<scalar>(os, "functions<scalar>", dict_, scalarFuncs_);
    writeFuncsImpl<vector>(os, "functions<vector>", dict_, vectorFuncs_);
}


// ************************************************************************* //

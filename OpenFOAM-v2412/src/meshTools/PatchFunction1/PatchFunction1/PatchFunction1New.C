/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "ConstantField.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::PatchFunction1<Type>>
Foam::PatchFunction1<Type>::New
(
    const polyPatch& pp,
    const word& entryName,
    const entry* eptr,
    const dictionary& dict,
    const bool faceValues,
    const bool mandatory
)
{
    word modelType;

    const dictionary* coeffs = (eptr ? eptr->dictPtr() : nullptr);

    if (coeffs)
    {
        // Dictionary entry

        DebugInFunction
            << "For " << entryName << " with dictionary entries: "
            << flatOutput(coeffs->toc()) << nl;

        // The "type" entry - mandatory, no fallback type provided
        coeffs->readEntry
        (
            "type",
            modelType,
            keyType::LITERAL,
            IOobjectOption::MUST_READ
        );
    }
    else if (eptr)
    {
        // Primitive entry
        // - word : the modelType, or uniform/nonuniform
        // - non-word : value for constant (uniform) function

        DebugInFunction
            << "For " << entryName << " with primitive entry" << nl;

        ITstream& is = eptr->stream();

        if (is.peek().isWord())
        {
            modelType = is.peek().wordToken();
        }
        else
        {
            // A value - compatibility for reading uniform (constant) field

            const Type constValue = pTraits<Type>(is);

            return autoPtr<PatchFunction1<Type>>
            (
                new PatchFunction1Types::ConstantField<Type>
                (
                    pp,
                    entryName,
                    constValue,
                    dict,
                    faceValues
                )
            );
        }

        // Looks like a normal field entry?
        if (modelType == "uniform" || modelType == "nonuniform")
        {
            return autoPtr<PatchFunction1<Type>>
            (
                new PatchFunction1Types::ConstantField<Type>
                (
                    pp,
                    eptr,
                    entryName,
                    dict,
                    faceValues
                )
            );
        }

        // Fallthrough
    }


    if (modelType.empty())
    {
        // Entry missing

        if (mandatory)
        {
            FatalIOErrorInFunction(dict)
                << "Missing or invalid PatchFunction1 entry: "
                << entryName << nl
                << exit(FatalIOError);
        }

        return nullptr;
    }
    else if (!coeffs)
    {
        // Primitive entry. Coeffs dictionary is optional.
        // Use keyword() - not entryName - for compatibility lookup!

        const word& kw =
        (
            eptr
          ? eptr->keyword()  // Could be a compatibility lookup
          : entryName
        );

        coeffs = &dict.optionalSubDict(kw + "Coeffs", keyType::LITERAL);
    }


    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown PatchFunction1 type "
            << modelType << " for " << entryName
            << "\n\nValid PatchFunction1 types :\n"
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return ctorPtr(pp, modelType, entryName, *coeffs, faceValues);
}


template<class Type>
Foam::autoPtr<Foam::PatchFunction1<Type>>
Foam::PatchFunction1<Type>::New
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues,
    const bool mandatory
)
{
    return PatchFunction1<Type>::New
    (
        pp,
        entryName,
        dict.findEntry(entryName, keyType::LITERAL),
        dict,
        faceValues,
        mandatory
    );
}


template<class Type>
Foam::autoPtr<Foam::PatchFunction1<Type>>
Foam::PatchFunction1<Type>::NewCompat
(
    const polyPatch& pp,
    const word& entryName,
    std::initializer_list<std::pair<const char*,int>> compat,
    const dictionary& dict,
    const bool faceValues,
    const bool mandatory
)
{
    return PatchFunction1<Type>::New
    (
        pp,
        entryName,
        dict.findCompat(entryName, compat, keyType::LITERAL),
        dict,
        faceValues,
        mandatory
    );
}


template<class Type>
Foam::autoPtr<Foam::PatchFunction1<Type>>
Foam::PatchFunction1<Type>::NewIfPresent
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
{
    // mandatory = false
    return PatchFunction1<Type>::New(pp, entryName, dict, faceValues, false);
}


template<class Type>
Foam::refPtr<Foam::PatchFunction1<Type>>
Foam::PatchFunction1<Type>::New
(
    HashPtrTable<PatchFunction1<Type>>& cache,

    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    enum keyType::option matchOpt,
    const bool faceValues,
    const bool mandatory
)
{
    // See corresponding comments in Function1::New (caching version)

    refPtr<PatchFunction1<Type>> fref;  // return value

    // Try for direct cache hit
    fref.cref(cache.get(entryName));

    if (fref)
    {
        return fref;
    }


    // Lookup from dictionary
    const entry* eptr = dict.findEntry(entryName, matchOpt);

    if (eptr)
    {
        // Use keyword (potentially a wildcard) instead of entry name
        const auto& kw = eptr->keyword();

        // Try for a cache hit
        fref.cref(cache.get(kw));

        if (!fref)
        {
            // Create new entry
            auto fauto
            (
                PatchFunction1<Type>::New
                (
                    pp,
                    kw,
                    eptr,  // Already resolved
                    dict,
                    faceValues,
                    mandatory
                )
            );

            if (fauto)
            {
                // Cache the newly created function
                fref.cref(fauto.get());
                cache.set(kw, fauto);
            }
        }
    }

    if (mandatory && !fref)
    {
        FatalIOErrorInFunction(dict)
            << "No match for " << entryName << nl
            << exit(FatalIOError);
    }

    return fref;
}


// ************************************************************************* //

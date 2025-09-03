/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    const word& entryName,
    const entry* eptr,
    const dictionary& dict,
    const word& redirectType,
    const objectRegistry* obrPtr,
    const bool mandatory
)
{
    word modelType(redirectType);

    const dictionary* coeffs = (eptr ? eptr->dictPtr() : nullptr);

    if (coeffs)
    {
        // Dictionary entry

        DebugInFunction
            << "For " << entryName << " with dictionary entries: "
            << flatOutput(coeffs->toc()) << nl;

        // The "type" entry - mandatory if no redirectType provided
        coeffs->readEntry
        (
            "type",
            modelType,
            keyType::LITERAL,
            (
                modelType.empty()
              ? IOobjectOption::MUST_READ : IOobjectOption::READ_IF_PRESENT
            )
        );

        // Fallthrough
    }
    else if (eptr)
    {
        // Primitive entry
        // - word : the modelType
        // - non-word : value for constant function

        DebugInFunction
            << "For " << entryName << " with primitive entry" << nl;

        ITstream& is = eptr->stream();

        if (is.peek().isWord())
        {
            modelType = is.peek().wordToken();
        }
        else
        {
            // A value - compatibility for reading constant

            const Type constValue = pTraits<Type>(is);

            return autoPtr<Function1<Type>>
            (
                new Function1Types::Constant<Type>
                (
                    entryName,
                    constValue,
                    obrPtr
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
                << "Missing or invalid Function1 entry: "
                << entryName << nl
                << exit(FatalIOError);
        }

        return nullptr;
    }
    else if (!coeffs)
    {
        // Primitive entry. Coeffs dictionary is optional.

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
            << "Unknown Function1 type "
            << modelType << " for " << entryName
            << "\n\nValid Function1 types :\n"
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return ctorPtr(entryName, *coeffs, obrPtr);
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    const word& entryName,
    const dictionary& dict,
    const word& redirectType,
    const objectRegistry* obrPtr,
    const bool mandatory
)
{
    return Function1<Type>::New
    (
        entryName,
        dict.findEntry(entryName, keyType::LITERAL),
        dict,
        redirectType,
        obrPtr,
        mandatory
    );
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::NewCompat
(
    const word& entryName,
    std::initializer_list<std::pair<const char*,int>> compat,
    const dictionary& dict,
    const word& redirectType,
    const objectRegistry* obrPtr,
    const bool mandatory
)
{
    return Function1<Type>::New
    (
        entryName,
        dict.findCompat(entryName, compat, keyType::LITERAL),
        dict,
        redirectType,
        obrPtr,
        mandatory
    );
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr,
    const bool mandatory
)
{
    return Function1<Type>::New(entryName, dict, word::null, obrPtr, mandatory);
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::NewIfPresent
(
    const word& entryName,
    const dictionary& dict,
    const word& redirectType,
    const objectRegistry* obrPtr
)
{
    // mandatory = false
    return Function1<Type>::New(entryName, dict, redirectType, obrPtr, false);
}


template<class Type>
Foam::autoPtr<Foam::Function1<Type>>
Foam::Function1<Type>::NewIfPresent
(
    const word& entryName,
    const dictionary& dict,
    const objectRegistry* obrPtr
)
{
    // mandatory = false
    return Function1<Type>::New(entryName, dict, word::null, obrPtr, false);
}


template<class Type>
Foam::refPtr<Foam::Function1<Type>>
Foam::Function1<Type>::New
(
    HashPtrTable<Function1<Type>>& cache,

    const word& entryName,
    const dictionary& dict,
    enum keyType::option matchOpt,
    const objectRegistry* obrPtr,
    const bool mandatory
)
{
    // Use the dictionary to find the keyword (allowing wildcards).
    // Alternative would be to have
    // a HashTable where the key type uses a wildcard match


    refPtr<Function1<Type>> fref;  // return value

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
                Function1<Type>::New
                (
                    kw,
                    eptr,  // Already resolved
                    dict,
                    word::null,
                    obrPtr,
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


/// template<class Type>
/// Foam::refPtr<Foam::Function1<Type>>
/// Foam::Function1<Type>::NewOrDefault
/// (
///     HashPtrTable<Function1<Type>>& cache,
///
///     const word& entryName,
///     const dictionary& dict,
///     const Type& deflt,
///     enum keyType::option matchOpt,
///     const objectRegistry* obrPtr
/// )
/// {
///     auto fref
///     (
///         Function1<Type>::New(entryName, dict, cache, matchOpt, obrPtr,false)
///     );
///
///     if (!fref)
///     {
///         fref.reset(new Function1Types::Constant<Type>("default", deflt));
///     }
///
///     return fref;
/// }

// ************************************************************************* //

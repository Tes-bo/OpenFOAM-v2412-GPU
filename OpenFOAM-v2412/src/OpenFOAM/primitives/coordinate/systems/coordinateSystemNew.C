/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "objectRegistry.H"
#include "cartesianCS.H"
#include "indirectCS.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::New
(
    const word& modelType,
    const dictionary& dict,
    IOobjectOption::readOption readOrigin,
    const objectRegistry* obrPtr
)
{
    // Direct dispatch
    // - treat missing modelType as 'cartesian'

    if (modelType.empty())
    {
        return autoPtr<coordinateSystem>
        (
            new coordSystem::cartesian(dict, readOrigin)
        );
    }


    // Dispatch with objectRegistry reference (if possible)
    if (obrPtr)
    {
        auto* ctorPtr = registryConstructorTable(modelType);

        if (ctorPtr)
        {
            return autoPtr<coordinateSystem>
            (
                ctorPtr(*obrPtr, dict, readOrigin)
            );
        }
    }

    // Regular dispatch
    // Note: everything with a registry constructor also has a
    // dictionary constructor, so just need to print those on error.

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "coordinate system",
             modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<coordinateSystem>(ctorPtr(dict, readOrigin));
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::New
(
    const dictionary& dict,
    const word& dictName,
    IOobjectOption::readOption readOrigin,
    const objectRegistry* obrPtr
)
{
    const dictionary* dictPtr = nullptr;

    // If dictName is non-empty: treat as mandatory
    // Include fallback handling of 'coordinateSystem' sub-dictionary
    // In 1806 and earlier, this was handled (rather poorly) in the
    // coordinateSystem constructor itself

    if (!dictName.empty())
    {
        const auto finder = dict.csearch(dictName, keyType::LITERAL);

        if (finder.isDict())
        {
            dictPtr = finder.dictPtr();
        }
        else
        {
            // Missing or primitive entry: trigger fatal error
            dictPtr = &(dict.subDict(dictName, keyType::LITERAL));
        }
    }
    else
    {
        // Search for "coordinateSystem" sub-dictionary

        const auto finder =
            dict.csearch(coordinateSystem::typeName, keyType::LITERAL);

        if (finder.isDict())
        {
            dictPtr = finder.dictPtr();
        }
        else if (finder.isStream())
        {
            const word csName(finder.stream());

            // Deprecated, unsupported syntax
            if (error::master())
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Ignoring '" << coordinateSystem::typeName
                    << "' as a keyword. Perhaps you meant this instead?" << nl
                    << '{' << nl
                    << "    type " << coordSystem::indirect::typeName
                    << ';' << nl
                    << "    name " << csName << ';' << nl
                    << '}' << nl
                    << std::endl;

                error::warnAboutAge("syntax change", 1806);
            }
        }
    }


    if (dictPtr)
    {
        // Using a sub-dictionary
        // - the 'origin' can be optional
        readOrigin = IOobjectOption::lazierRead(readOrigin);
    }
    else
    {
        // Using top-level dictionary
        dictPtr = &dict;
    }


    // The coordinate-system type (if not cartesian)
    word modelType;
    dictPtr->readIfPresent("type", modelType, keyType::LITERAL);

    return coordinateSystem::New
    (
        modelType,
        *dictPtr,
        readOrigin,
        obrPtr
    );
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::NewIfPresent
(
    const dictionary& dict,
    const word& dictName,
    const objectRegistry* obrPtr
)
{
    const dictionary* dictPtr = nullptr;

    if
    (
        dictName.empty()
     || (dictPtr = dict.findDict(dictName, keyType::LITERAL)) == nullptr
    )
    {
        return nullptr;
    }

    // The coordinate-system type (if not cartesian)
    word modelType;
    dictPtr->readIfPresent("type", modelType, keyType::LITERAL);

    return coordinateSystem::New
    (
        modelType,
        *dictPtr,
        IOobjectOption::READ_IF_PRESENT,
        obrPtr
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::New
(
    Istream& is,
    IOobjectOption::readOption readOrigin
)
{
    const word csName(is);
    const dictionary dict(is);

    // The coordinate-system type (if not cartesian)
    word modelType;
    dict.readIfPresent("type", modelType, keyType::LITERAL);

    auto cs = coordinateSystem::New(modelType, dict, readOrigin);
    cs->rename(csName);

    return cs;
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::New
(
    const word& modelType,
    const objectRegistry& obr,
    const dictionary& dict,
    IOobjectOption::readOption readOrigin
)
{
    return New(modelType, dict, readOrigin, &obr);
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::New
(
    const word& modelType,
    const dictionary& dict,
    IOobjectOption::readOption readOrigin
)
{
    return New(modelType, dict, readOrigin, nullptr);
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& dictName,
    IOobjectOption::readOption readOrigin
)
{
    return New(dict, dictName, readOrigin, &obr);
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::New
(
    const dictionary& dict,
    const word& dictName,
    IOobjectOption::readOption readOrigin
)
{
    return New(dict, dictName, readOrigin, nullptr);
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::NewIfPresent
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& dictName
)
{
    return NewIfPresent(dict, dictName, &obr);
}


Foam::autoPtr<Foam::coordinateSystem>
Foam::coordinateSystem::NewIfPresent
(
    const dictionary& dict,
    const word& dictName
)
{
    return NewIfPresent(dict, dictName, nullptr);
}


// ************************************************************************* //

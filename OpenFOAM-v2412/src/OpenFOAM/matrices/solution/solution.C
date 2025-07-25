/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "solution.H"
#include "HashPtrTable.H"
#include "Function1.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineDebugSwitchWithName(solution, "solution", 0);
    registerDebugSwitchWithName(solution, solution, "solution");
}

// List of sub-dictionaries to rewrite
static const Foam::List<Foam::word> subDictNames
({
    "preconditioner",
    "smoother"
});


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solution::read(const dictionary& dict)
{
    const dictionary* dictptr;

    if ((dictptr = dict.findDict("cache")) != nullptr)
    {
        cache_ = *dictptr;
        caching_ = cache_.getOrDefault("active", true);
    }

    if ((dictptr = dict.findDict("relaxationFactors")) != nullptr)
    {
        const dictionary& relaxDict = *dictptr;

        bool needsCompat = true;

        if ((dictptr = relaxDict.findDict("fields")) != nullptr)
        {
            needsCompat = false;
            fieldRelaxDict_ = *dictptr;
            fieldRelaxCache_.clear();
        }

        if ((dictptr = relaxDict.findDict("equations")) != nullptr)
        {
            needsCompat = false;
            eqnRelaxDict_ = *dictptr;
            eqnRelaxCache_.clear();
        }

        if (needsCompat)
        {
            // backwards compatibility
            fieldRelaxDict_.clear();
            fieldRelaxCache_.clear();

            for (const word& e : relaxDict.toc())
            {
                scalar value = relaxDict.get<scalar>(e);

                if (e.starts_with('p'))
                {
                    fieldRelaxDict_.add(e, value);
                }
                else if (e.starts_with("rho"))
                {
                    fieldRelaxDict_.add(e, value);
                }
            }

            eqnRelaxDict_ = relaxDict;
            eqnRelaxCache_.clear();
        }


        fieldRelaxDefault_ = Function1<scalar>::NewIfPresent
        (
            "default",
            fieldRelaxDict_,
            &db()
        );
        if (!fieldRelaxDefault_)
        {
            fieldRelaxDefault_.reset
            (
                new Function1Types::Constant<scalar>("default", 0, &db())
            );
        }

        eqnRelaxDefault_ = Function1<scalar>::NewIfPresent
        (
            "default",
            eqnRelaxDict_,
            &db()
        );
        if (!eqnRelaxDefault_)
        {
            eqnRelaxDefault_.reset
            (
                new Function1Types::Constant<scalar>("default", 0, &db())
            );
        }

        DebugInfo
            << "Relaxation factors:" << nl
            << "fields: " << fieldRelaxDict_ << nl
            << "equations: " << eqnRelaxDict_ << endl;
    }

    if ((dictptr = dict.findDict("solvers")) != nullptr)
    {
        solvers_ = *dictptr;
        upgradeSolverDict(solvers_);
    }
}


const Foam::dictionary& Foam::solution::selectedDict() const
{
    word select;

    if (readIfPresent("select", select, keyType::LITERAL))
    {
        return subDict(select);
    }

    return *this;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solution::solution
(
    const objectRegistry& obr,
    IOobjectOption::readOption rOpt,
    const fileName& dictName,
    const dictionary* fallback
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            rOpt,
            IOobjectOption::NO_WRITE
        ),
        fallback
    ),
    cache_(),
    caching_(false),
    fieldRelaxDict_(),
    eqnRelaxDict_(),
    solvers_()
{
    // Treat as READ_MODIFIED whenever possible
    if
    (
        readOpt() == IOobjectOption::MUST_READ
     || (isReadOptional() && headerOk())
    )
    {
        readOpt(IOobjectOption::READ_MODIFIED);
        addWatch();
    }

    // Update: from values read or copied in
    if
    (
        readOpt() == IOobjectOption::READ_MODIFIED
     || !dictionary::empty()
    )
    {
        read(selectedDict());
    }
}


Foam::solution::solution
(
    const objectRegistry& obr,
    const fileName& dictName,
    const dictionary* fallback
)
:
    solution(obr, obr.readOpt(), dictName, fallback)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// No default destructor in header (incomplete types)
Foam::solution::~solution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::solution::upgradeSolverDict
(
    dictionary& dict,
    const bool verbose
)
{
    label nChanged = 0;

    // backward compatibility:
    // recast primitive entries into dictionary entries
    for (const entry& dEntry : dict)
    {
        if (dEntry.isStream())
        {
            ITstream& is = dEntry.stream();
            word name(is);
            dictionary subdict;

            subdict.add("solver", name);
            subdict <<= dictionary(is);

            // preconditioner and smoother entries can be
            // 1) primitiveEntry w/o settings,
            // 2) or a dictionaryEntry.
            // transform primitiveEntry with settings -> dictionaryEntry
            for (const word& dictName : subDictNames)
            {
                ITstream* streamPtr =
                    subdict.findStream(dictName, keyType::LITERAL);

                if (streamPtr)
                {
                    auto& is = *streamPtr;
                    is >> name;

                    if (!is.eof())
                    {
                        dictionary newDict;
                        newDict.add(dictName, name);
                        newDict <<= dictionary(is);

                        subdict.set(dictName, newDict);
                    }
                }
            }

            // write out information to help people adjust to the new syntax
            if (verbose && Pstream::master())
            {
                Info<< "// using new solver syntax:\n"
                    << dEntry.keyword() << subdict << endl;
            }

            // overwrite with dictionary entry
            dict.set(dEntry.keyword(), subdict);

            ++nChanged;
        }
    }

    return nChanged;
}


bool Foam::solution::cache(const word& name) const
{
    if (caching_)
    {
        DebugInfo<< "Cache: find entry for " << name << endl;
        return cache_.found(name);
    }

    return false;
}


// void Foam::solution::enableCache(const word& name) const
// {
//     if (!cache_.found(name))
//     {
//         DebugInfo<< "Cache: enable cache for " << name << endl;
//         cache_.add(name, true);
//         caching_ = true;
//     }
// }


bool Foam::solution::relaxField(const word& name) const
{
    DebugInfo
        << "Field relaxation factor for " << name
        << " is " << (fieldRelaxDict_.found(name) ? "set" : "unset") << endl;

    return fieldRelaxDict_.found(name) || fieldRelaxDict_.found("default");
}


bool Foam::solution::relaxEquation(const word& name) const
{
    DebugInfo<< "Find equation relaxation factor for " << name << endl;
    return eqnRelaxDict_.found(name) || eqnRelaxDict_.found("default");
}


bool Foam::solution::relaxField(const word& name, scalar& factor) const
{
    DebugInfo<< "Lookup field relaxation factor for " << name << endl;

    if (fieldRelaxDict_.found(name))
    {
        factor = Function1<scalar>::New
        (
            fieldRelaxCache_,  // cache
            name,
            fieldRelaxDict_,
            keyType::REGEX,
            &db()
        )().value(time().timeOutputValue());

        DebugInfo
            << "Field relaxation factor for " << name
            << " is " << factor
            << " (from Function1)" << endl;

        return true;
    }
    else if (fieldRelaxDict_.found("default") && fieldRelaxDefault_)
    {
        factor = fieldRelaxDefault_->value(time().timeOutputValue());

        DebugInfo
            << "Field relaxation factor for " << name
            << " is " << factor
            << " (from default " << eqnRelaxDefault_->type() << ')' << endl;

        return true;
    }

    // Fallthrough - nothing found

    DebugInfo<< "No field relaxation factor for " << name << endl;

    return false;
}


bool Foam::solution::relaxEquation(const word& name, scalar& factor) const
{
    DebugInfo<< "Lookup equation relaxation factor for " << name << endl;

    if (eqnRelaxDict_.found(name))
    {
        factor = Function1<scalar>::New
        (
            eqnRelaxCache_,  // cache
            name,
            eqnRelaxDict_,
            keyType::REGEX,
            &db()
        )().value(time().timeOutputValue());

        DebugInfo
            << "Equation relaxation factor for " << name
            << " is " << factor
            << " (from Function1)" << endl;

        return true;
    }
    else if (eqnRelaxDict_.found("default") && eqnRelaxDefault_)
    {
        factor = eqnRelaxDefault_->value(time().timeOutputValue());

        DebugInfo
            << "Equation relaxation factor for " << name
            << " is " << factor
            << " (from default " << eqnRelaxDefault_->type() << ')' << endl;

        return true;
    }

    // Fallthrough - nothing found

    DebugInfo<< "No equation relaxation factor for " << name << endl;

    return false;
}


Foam::scalar Foam::solution::fieldRelaxationFactor(const word& name) const
{
    // Any initial value
    scalar factor = 0;

    if (!relaxField(name, factor))
    {
        FatalIOErrorInFunction(fieldRelaxDict_)
            << "Cannot find field relaxation factor for '" << name
            << "' or a suitable default value." << nl
            << exit(FatalIOError);
    }

    return factor;
}


Foam::scalar Foam::solution::equationRelaxationFactor(const word& name) const
{
    // Any initial value
    scalar factor = 0;

    if (!relaxEquation(name, factor))
    {
        FatalIOErrorInFunction(eqnRelaxDict_)
            << "Cannot find equation relaxation factor for '" << name
            << "' or a suitable default value."
            << exit(FatalIOError);
    }

    return factor;
}


const Foam::dictionary& Foam::solution::solutionDict() const
{
    return selectedDict();
}


const Foam::dictionary& Foam::solution::solutionDict(const word& name) const
{
    DebugInfo<< "Lookup subDict : " << name << endl;
    return selectedDict().subDict(name);
}


const Foam::dictionary& Foam::solution::solversDict() const
{
    return solvers_;
}


const Foam::dictionary& Foam::solution::solverDict(const word& name) const
{
    DebugInfo<< "Lookup solver for " << name << endl;
    return solvers_.subDict(name);
}


const Foam::dictionary& Foam::solution::solver(const word& name) const
{
    DebugInfo<< "Lookup solver for " << name << endl;
    return solvers_.subDict(name);
}


bool Foam::solution::read()
{
    if (regIOobject::read())
    {
        read(selectedDict());

        return true;
    }

    return false;
}


// ************************************************************************* //

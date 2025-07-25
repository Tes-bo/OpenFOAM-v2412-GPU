/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
#include "expressionEntry.H"
#include "stringOps.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

    defineTypeNameAndDebug(exprDriver, 0);

} // End namespace expressions
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

int Foam::expressions::exprDriver::getSearchControls(const dictionary& dict)
{
    int val = 0;

    if (dict.getOrDefault("searchInMemory", true))
    {
        val |= int(searchControls::SEARCH_REGISTRY);
    }
    if (dict.getOrDefault("searchFiles", false))
    {
        val |= int(searchControls::SEARCH_FILES);
    }
    if (dict.getOrDefault("cacheReadFields", false))
    {
        val |= int(searchControls::CACHE_READ_FIELDS);
    }

    return val;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
#if 0
static string getEntryString
(
    const dictionary& dict,
    const string& key
)
{
    const entry* eptr = dict.findEntry(key, keyType::REGEX_RECURSIVE);

    if (!eptr)
    {
        FatalErrorInFunction
            << "Entry " << key << " not found in "
            << dict.name() << nl
            << exit(FatalError);
    }
    else if (eptr->isDict())
    {
        FatalErrorInFunction
            << "Entry " << key << " found in "
            << dict.name() << " but is a dictionary" << nl
            << exit(FatalError);
    }

    return exprTools::expressionEntry::evaluate(*eptr);
}
#endif


template<class Type>
static void shallowCloneFunctions
(
    HashTable<refPtr<Function1<Type>>>& dest,
    const HashTable<refPtr<Function1<Type>>>& rhs
)
{
    // Add in shallow copy for other functions
    forAllConstIters(rhs, iter)
    {
        const word& key = iter.key();

        if (!dest.found(key))
        {
            refPtr<Function1<Type>> func;
            func.cref(iter.val().shallowClone());

            dest.emplace_set(key, std::move(func));
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::expressions::exprDriver::resetTimeReference(const TimeState* ts)
{
    timeStatePtr_ = ts;
}


void Foam::expressions::exprDriver::resetTimeReference(const TimeState& ts)
{
    timeStatePtr_ = &ts;
}


void Foam::expressions::exprDriver::resetDb(const objectRegistry* obrPtr)
{
    obrPtr_ = obrPtr;

    forAllIters(scalarFuncs_, iter)
    {
        auto& funcPtr = iter.val();
        if (funcPtr && !funcPtr.is_const())
        {
            (*funcPtr).resetDb(obrPtr_);
        }
    }
    forAllIters(vectorFuncs_, iter)
    {
        auto& funcPtr = iter.val();
        if (funcPtr && !funcPtr.is_const())
        {
            (*funcPtr).resetDb(obrPtr_);
        }
    }
}


void Foam::expressions::exprDriver::resetDb(const objectRegistry& db)
{
    resetDb(&db);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::exprDriver::exprDriver
(
    enum searchControls search,
    const dictionary& dict
)
:
    dict_(dict),
    result_(),
    variableStrings_(),
    variables_(16),
    scalarFuncs_(0),
    vectorFuncs_(0),
    contextObjects_(0),
    arg1Value_(0),
    timeStatePtr_(nullptr),
    obrPtr_(nullptr),
    stashedTokenId_(0),

    // Controls
    debugScanner_(dict.getOrDefault("debug.scanner", false)),
    debugParser_(dict.getOrDefault("debug.parser", false)),
    allowShadowing_(dict.getOrDefault("allowShadowing", false)),
    prevIterIsOldTime_(dict.getOrDefault("prevIterIsOldTime", false)),
    searchCtrl_(search)
{}


Foam::expressions::exprDriver::exprDriver
(
    const exprDriver& rhs,
    const dictionary& dict
)
:
    dict_(dict),
    result_(rhs.result_),
    variableStrings_(rhs.variableStrings_),
    variables_(rhs.variables_),
    scalarFuncs_(0),
    vectorFuncs_(0),
    contextObjects_(rhs.contextObjects_),
    arg1Value_(rhs.arg1Value_),
    timeStatePtr_(rhs.timeStatePtr_),
    obrPtr_(rhs.obrPtr_),
    stashedTokenId_(0),

    // Controls
    debugScanner_(rhs.debugScanner_),
    debugParser_(rhs.debugParser_),
    allowShadowing_(rhs.allowShadowing_),
    prevIterIsOldTime_(rhs.prevIterIsOldTime_),

    searchCtrl_(rhs.searchCtrl_)
{
    // Partially like readDict()

    // Create Function1s from dictionary content
    resetFunctions(dict_);

    // Add in shallow copy for other functions
    shallowCloneFunctions(scalarFuncs_, rhs.scalarFuncs_);
    shallowCloneFunctions(vectorFuncs_, rhs.vectorFuncs_);
}


Foam::expressions::exprDriver::exprDriver
(
    const dictionary& dict
)
:
    exprDriver
    (
        searchControls(exprDriver::getSearchControls(dict)),
        dict
    )
{
    readDict(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::TimeState* Foam::expressions::exprDriver::timeState() const noexcept
{
    if (timeStatePtr_)
    {
        return timeStatePtr_;
    }
    else if (obrPtr_)
    {
        return &(obrPtr_->time());
    }
    return nullptr;
}


Foam::scalar Foam::expressions::exprDriver::timeValue() const
{
    if (timeStatePtr_)
    {
        return timeStatePtr_->value();
    }
    else if (obrPtr_)
    {
        return obrPtr_->time().value();
    }
    return 0;
}


Foam::scalar Foam::expressions::exprDriver::deltaT() const
{
    if (timeStatePtr_)
    {
        return timeStatePtr_->deltaTValue();
    }
    else if (obrPtr_)
    {
        return obrPtr_->time().deltaTValue();
    }
    return 0;
}


bool Foam::expressions::exprDriver::readDict
(
    const dictionary& dict
)
{
    dict.readIfPresent("debug.driver", debug);

    // Regular variables (optional)
    variableStrings_ = readVariableStrings(dict, "variables", false);

    // Create Function1s from dictionary content
    resetFunctions(dict);

    // readTable("lookuptables2D", dict, lookup2D_);

    return true;
}


void Foam::expressions::exprDriver::clearResult()
{
    result_.clear();
}


bool Foam::expressions::exprDriver::update()
{
    return true;
}


void Foam::expressions::exprDriver::updateSpecialVariables(bool force)
{}


void Foam::expressions::exprDriver::clearVariables()
{
    variables_.clear();
    addVariables(variableStrings_, false);
}


void Foam::expressions::exprDriver::evaluateVariable
(
    const word& varName,
    const expressions::exprString& expr
)
{
    parse(expr);
    result_.testIfSingleValue();

    DebugInfo
        << "Evaluating: " << expr << " -> " << varName << endl
        << result_;

    // Overwrite with a copy
    variables_.set(varName, exprResult(result_));
}


void Foam::expressions::exprDriver::evaluateVariableRemote
(
    string remote,
    const word& varName,
    const expressions::exprString& expr
)
{
    NotImplemented;
}


Foam::expressions::exprResult
Foam::expressions::exprDriver::getRemoteResult
(
    const exprDriver& other
) const
{
    // With warnings (noWarn = false)
    return other.result().getUniform(this->size(), false);
}


void Foam::expressions::exprDriver::addVariables
(
    const expressions::exprString& expr,
    bool clear
)
{
    if (clear)
    {
        clearVariables();
    }

    // Allow inline list of semicolon-separated variables
    const auto varExpressions =
        stringOps::split<expressions::exprString>(expr, ';');

    for (const auto& subMatch : varExpressions)
    {
        string varExpr(stringOps::trim(subMatch.str()));
        if (varExpr.empty())
        {
            continue;
        }

        // Split on '=' for lhsExpr = rhsExpr
        //
        // varName = rhsExpr
        // varName{where} = rhsExpr

        const auto eqPos = varExpr.find('=');

        if (eqPos == std::string::npos)
        {
            FatalIOErrorInFunction(dict_)
                << "No '=' found in expression " << varExpr << nl << nl
                << exit(FatalIOError);
        }

        // The RHS
        expressions::exprString rhsExpr
        (
            expressions::exprString::toExpr
            (
                stringOps::trim(varExpr.substr(eqPos+1))
            )
        );

        // The LHS
        varExpr.resize(eqPos);
        stringOps::inplaceTrim(varExpr);

        // Check for varName{where}
        const auto lbrace = varExpr.find('{');

        if (lbrace != std::string::npos)
        {
            const auto rbrace = varExpr.find('}');

            if (rbrace == std::string::npos || rbrace < lbrace)
            {
                FatalErrorInFunction
                    // << "Context: " << driverContext_ << nl
                    << "No closing '}' found in " << varExpr << nl
                    << exit(FatalError);
            }
            else if (lbrace+1 == rbrace)
            {
                FatalErrorInFunction
                    // << "Context: " << driverContext_ << nl
                    << "Empty '{}' location in " << varExpr << nl
                    << exit(FatalError);
            }

            const word varName(word::validate(varExpr.substr(0, lbrace)));

            const expressions::exprString remoteExpr
            (
                expressions::exprString::toExpr
                (
                    varExpr.substr(lbrace+1, rbrace-lbrace-1)
                )
            );

            // Fails if derived class does not implement!

            evaluateVariableRemote(remoteExpr, varName, rhsExpr);
        }
        else
        {
            const word varName(word::validate(varExpr));

            evaluateVariable(varName, rhsExpr);
        }
    }
}


void Foam::expressions::exprDriver::addVariables
(
    const UList<expressions::exprString>& list,
    bool clear
)
{
    if (clear)
    {
        clearVariables();
    }

    for (const auto& expr : list)
    {
        addVariables(expr, false); // No clear (already done)
    }
}


void Foam::expressions::exprDriver::setDebugging
(
    bool scannerDebug,
    bool parserDebug
)
{
    debugScanner_ = scannerDebug;
    debugParser_ = parserDebug;
}


void Foam::expressions::exprDriver::setDebugging
(
    const exprDriver& rhs
)
{
    debugScanner_ = rhs.debugScanner_;
    debugParser_ = rhs.debugParser_;
}


bool Foam::expressions::exprDriver::setCaching(bool on) noexcept
{
    int val(searchCtrl_);
    bool old = (val & searchControls::CACHE_READ_FIELDS);

    if (!on)
    {
        // Off
        val &= ~(searchControls::CACHE_READ_FIELDS);
    }
    else if (!old)
    {
        // Toggled on.
        // Caching read fields implies both registry and disk use
        val |=
        (
            searchControls::SEARCH_REGISTRY
          | searchControls::SEARCH_FILES
          | searchControls::CACHE_READ_FIELDS
        );
    }

    searchCtrl_ = searchControls(val);

    return old;
}


void Foam::expressions::exprDriver::setSearchBehaviour
(
    enum searchControls search,
    const bool caching
)
{
    int val(search);
    if (caching || (val & searchControls::CACHE_READ_FIELDS))
    {
        // Caching read fields implies both registry and disk use
        val |=
        (
            searchControls::SEARCH_REGISTRY
          | searchControls::SEARCH_FILES
          | searchControls::CACHE_READ_FIELDS
        );
    }
    searchCtrl_ = searchControls(val);

    #ifdef FULLDEBUG
    Info<< "Searching "
        << " registry:" << searchRegistry()
        << " disk:" << searchFiles()
        << " cache-read:" << cacheReadFields() << nl;
    #endif
}


void Foam::expressions::exprDriver::setSearchBehaviour
(
    const exprDriver& rhs
)
{
    searchCtrl_ = rhs.searchCtrl_;
}


// ************************************************************************* //

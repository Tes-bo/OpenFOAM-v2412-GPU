/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "fvExprDriver.H"
#include "fvExprDriverWriter.H"
#include "expressionEntry.H"
#include "exprResultGlobals.H"

#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

    defineTypeNameAndDebug(fvExprDriver, 0);
    defineRunTimeSelectionTable(fvExprDriver, dictionary);
    defineRunTimeSelectionTable(fvExprDriver, idName);

} // End namespace expressions
} // End namespace Foam

// Currently not working?
bool Foam::expressions::fvExprDriver::cacheSets_ = true;

const Foam::fvMesh* Foam::expressions::fvExprDriver::defaultMeshPtr_ = nullptr;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::fvMesh& Foam::expressions::fvExprDriver::defaultMesh()
{
    if (!defaultMeshPtr_)
    {
        FatalErrorInFunction
            << "No default mesh set" << nl
            << "Try the 'fvExprDriverFunctionObject' as a workaround"
            << endl
            << abort(FatalError);
    }

    return *defaultMeshPtr_;
}


const Foam::fvMesh* Foam::expressions::fvExprDriver::resetDefaultMesh
(
    const fvMesh& mesh,
    const bool force
)
{
    const fvMesh* ptr = defaultMeshPtr_;

    if (force || (ptr != nullptr))
    {
        defaultMeshPtr_ = &mesh;
    }

    return ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::expressions::fvExprDriver::fvExprDriver
(
    enum exprDriver::searchControls search,
    const dictionary& dict
)
:
    expressions::exprDriver(search, dict),
    globalScopes_(),
    delayedVariables_(),
    storedVariables_(),
    specialVariablesIndex_(-1),
    otherMeshName_(),
    writer_(nullptr)
{}


Foam::expressions::fvExprDriver::fvExprDriver
(
    const fvExprDriver& rhs,
    const dictionary& dict
)
:
    expressions::exprDriver(rhs, dict),
    globalScopes_(rhs.globalScopes_),
    delayedVariables_(rhs.delayedVariables_),
    storedVariables_(rhs.storedVariables_),
    specialVariablesIndex_(rhs.specialVariablesIndex_),
    otherMeshName_(),
    writer_(nullptr)
{}


Foam::expressions::fvExprDriver::fvExprDriver
(
    const dictionary& dict
)
:
    expressions::exprDriver(dict),
    globalScopes_(),
    delayedVariables_(),
    storedVariables_(),
    specialVariablesIndex_(-1),
    otherMeshName_(),
    writer_(nullptr)
{
    readDict(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expressions::fvExprDriver::~fvExprDriver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::fvExprDriver::readDict
(
    const dictionary& dict
)
{
    expressions::exprDriver::readDict(dict);

    // fileNameList plugins;
    // if (dict.readIfPresent("functionPlugins", plugins))
    // {
    //     for (const fileName& libName : plugins)
    //     {
    //         this->mesh().time().libs().open
    //         (
    //             "libswak" + libName + "FunctionPlugin"  // verbose = true
    //         );
    //     }
    // }

    dict.readIfPresent("globalScopes", globalScopes_);

    const entry* eptr = nullptr;

    // Special variables

    if
    (
        // storedVariables
        (eptr = dict.findEntry("storedVariables", keyType::LITERAL))
     != nullptr
    )
    {
        ITstream& is = eptr->stream();

        if (writer_ && !storedVariables_.empty())
        {
            WarningInFunction
                // << "Context: " << driverContext_ << nl
                << "The 'storedVariables' was already read."
                << " No update from " << is
                << endl;
        }
        else
        {
            storedVariables_ = List<exprResultStored>(is);

            // Check for excess tokens
            dict.checkITstream(is, "storedVariables");
        }
    }

    if
    (
        // delayedVariables
        (eptr = dict.findEntry("delayedVariables", keyType::LITERAL))
     != nullptr
    )
    {
        ITstream& is = eptr->stream();

        if (writer_ && !delayedVariables_.empty())
        {
            WarningInFunction
                // << "Context: " << driverContext_ << nl
                << "Seems like 'delayedVariables' was already read."
                << " No update from " << is
                << endl;
        }
        else
        {
            List<exprResultDelayed> inputs(is);

            // Check for excess tokens
            dict.checkITstream(is, "delayedVariables");

            for (auto& var : inputs)
            {
                delayedVariables_.insert(var.name(), var);
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::fvExprDriver::updateSpecialVariables(bool force)
{
    const bool updated = this->update();

    const label eventIndex = mesh().time().timeIndex();
    const scalar eventTime = mesh().time().value();

    DebugInfo
        << "fvExprDriver::updateSpecialVariables(force="
        << force << ") Updated: " << updated << endl;

    if (specialVariablesIndex_ < 0)
    {
        DebugInfo
            << "First update: " << eventIndex << endl;

        specialVariablesIndex_ = eventIndex;

        for (exprResultStored& v : storedVariables_)
        {
            DebugInfo
                << v.name() << " = " << v.initialValueExpression()
                << " (has value "
                << v.hasValue() << ")" << endl;

            if (!v.hasValue())
            {
                DebugInfo
                    << "First value: " << v.initialValueExpression()
                    << " -> " << v.name() << endl;

                parse(v.initialValueExpression());
                v = result_;
                DebugInfo
                    << "Parser size: " << this->size() << nl
                    << "Calculated: " << result_ << nl
                    << "Stored: " << v << nl;
            }
        }
    }

    if (force || specialVariablesIndex_ != eventIndex)
    {
        DebugInfo
            << "Store variables: " << force << ' '
            << specialVariablesIndex_ << ' '
            << eventIndex << endl;

        for (exprResultStored& v : storedVariables_)
        {
            if (variables_.found(v.name()))
            {
                DebugInfo
                    << "Storing variable: " << v.name() << " "
                    << variables_[v.name()] << endl;

                v = variables_[v.name()];
            }
        }
        specialVariablesIndex_ = eventIndex;
    }

    forAllIters(delayedVariables_, iter)
    {
        DebugInfo
            << "Updating delayed variable " << iter().name() << endl;

        if (!iter().updateReadValue(eventTime))
        {
            const exprString& expr = iter().startupValueExpression();

            DebugInfo
                << "Evaluate: " << expr << endl;

            parse(expr);
            iter().setReadValue(result_);

            DebugInfo
                << "Value " << iter() << nl
                << "Type " << iter().valueType() << "("
                << result_.valueType() << ")" << endl;
        }
        else
        {
            DebugInfo
                << iter().name() << " updated without problem" << endl;
        }
    }
}


void Foam::expressions::fvExprDriver::clearVariables()
{
    DebugInfo
        << "Clearing variables" << endl;

    const scalar eventTime = mesh().time().value();

    (void)this->update();

    updateSpecialVariables();
    variables_.clear();
    for (exprResultStored& v : storedVariables_)
    {
        variables_.insert(v.name(), v);
    }

    addVariables(variableStrings_, false);

    forAllIters(delayedVariables_, iter)
    {
        iter().storeValue(eventTime);
    }
}


void Foam::expressions::fvExprDriver::evaluateVariable
(
    const word& varName,
    const expressions::exprString& expr
)
{
    const regIOobject* objPtr = mesh().findObject<regIOobject>(varName);

    if (!allowShadowing_ && objPtr)
    {
        WarningInFunction
            // << "Context: " << driverContext_ << nl
            << "Field '" << varName << "' (type " << objPtr->headerClassName()
            << ") is shadowed by a variable of the same name." << nl
            << "This may lead to trouble" << nl
            << "If this is OK set 'allowShadowing'"
            << " in the relevant parser" << nl
            << endl;
    }

    parse(expr);
    result_.testIfSingleValue();

    DebugInfo
        << "Evaluating: " << expr << " -> " << varName << endl
        << result_;


    // Assign
    if (delayedVariables_.found(varName))
    {
        // Avoid potential conflicts?
        variables_.erase(varName);

        DebugInfo
            << varName << " is delayed" << endl;

        // Copy assignment
        delayedVariables_[varName] = result_;
    }
    else
    {
        // Overwrite with a copy
        variables_.set(varName, exprResult(result_));
    }
}


void Foam::expressions::fvExprDriver::evaluateVariableRemote
(
    string remote,
    const word& varName,
    const expressions::exprString& expr
)
{
    DebugInfo
        << "Evaluating remote " << remote.c_str()
        << " : " << expr << " -> " << varName << endl;

    word driverType("patch");  // default is patch
    word identName, regionName;

    const auto slashPos = remote.find('/');
    if (slashPos != std::string::npos)
    {
        regionName = word::validate(remote.substr(slashPos+1));
        remote.resize(slashPos);
    }

    const auto quotePos = remote.find('\'');
    if (quotePos != std::string::npos)
    {
        driverType = word::validate(remote.substr(0, quotePos));
        identName = word::validate(remote.substr(quotePos+1));
    }
    else
    {
        identName = word::validate(remote);
    }

    if
    (
        driverType == "patch"
     &&
        (
            identName.empty()
         || identName == "volume"
         || identName == "internalField"
        )
    )
    {
        driverType = "internalField";
    }

    const fvMesh* pRegion = &(this->mesh());

    if (!regionName.empty())
    {
        pRegion = pRegion->time().cfindObject<fvMesh>(regionName);

        if (!pRegion)
        {
            FatalErrorInFunction
                << "Cannot resolve mesh region: " << regionName << nl
                << exit(FatalError);
        }
    }

    DebugInfo
        << "Call other with ("
        << driverType << ", " << identName << ", " << regionName << ")\n";

    autoPtr<fvExprDriver> otherDriver =
        fvExprDriver::New(driverType, identName, *pRegion);

    otherDriver->setSearchBehaviour(*this);
    otherDriver->setGlobalScopes(this->globalScopes_);

    otherDriver->parse(expr);

    exprResult otherResult(this->getRemoteResult(*otherDriver));

    // Check / re-check for uniform. Not normally needed
    if (!otherResult.isUniform())
    {
        otherResult.testIfSingleValue();
    }

    DebugInfo
        << "Remote result: " << otherResult << nl;

    // Assign
    if (delayedVariables_.found(varName))
    {
        // Avoid potential conflicts?
        variables_.erase(varName);

        DebugInfo
            << varName << " is delayed - setting" << nl;

        // Move assignment
        delayedVariables_[varName] = std::move(otherResult);
    }
    else
    {
        // Overwrite with a copy
        variables_.set(varName, std::move(otherResult));
    }
}


const Foam::fvMesh&
Foam::expressions::fvExprDriver::regionMesh
(
    const dictionary& dict,
    const fvMesh& mesh,
    bool readIfNecessary
)
{
    word regionName;

    if (!dict.readIfPresent("region", regionName))
    {
        DebugInFunction << "Using original mesh " << nl;
        return mesh;
    }

    DebugInFunction << "Using mesh " << regionName  << endl;

    fvMesh* meshPtr = mesh.time().getObjectPtr<fvMesh>(regionName);

    if (!meshPtr && readIfNecessary)
    {
        WarningInFunction
            << "Region " << regionName
            << " not in memory. Loading it" << endl;

        meshPtr = new fvMesh
        (
            IOobject
            (
                regionName,
                mesh.time().constant(),
                mesh.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::REGISTER
            )
        );

        meshPtr->polyMesh::store();
    }

    if (!meshPtr)
    {
        FatalErrorInFunction
            << "No mesh region loaded: " << regionName
            << endl;
    }

    return *meshPtr;
}


Foam::word Foam::expressions::fvExprDriver::getTypeOfField
(
    const word& fieldName
) const
{
    return getHeaderClassName(this->mesh(), fieldName);
}


Foam::word Foam::expressions::fvExprDriver::getFieldClassName
(
    const word& name
) const
{
    if (searchRegistry())
    {
        const regIOobject* ioptr = this->mesh().findObject<regIOobject>(name);

        if (ioptr)
        {
            return ioptr->type();
        }
    }

    if (searchFiles())
    {
        return getHeaderClassName(this->mesh(), name);
    }

    return word::null;
}


Foam::topoSetSource::sourceType
Foam::expressions::fvExprDriver::topoSetType(const word& setName) const
{
    IOobject io(topoSet::findIOobject(mesh(), setName));

    if (io.isHeaderClass<cellSet>())
    {
        return topoSetSource::sourceType::CELLSET_SOURCE;
    }
    if (io.isHeaderClass<faceSet>())
    {
        return topoSetSource::sourceType::FACESET_SOURCE;
    }
    if (io.isHeaderClass<pointSet>())
    {
        return topoSetSource::sourceType::POINTSET_SOURCE;
    }

    return topoSetSource::sourceType::UNKNOWN_SOURCE;
}


Foam::topoSetSource::sourceType
Foam::expressions::fvExprDriver::topoZoneType(const word& setName) const
{
    if (mesh().cellZones().findZoneID(setName) >= 0)
    {
        return topoSetSource::sourceType::CELLZONE_SOURCE;
    }

    if (mesh().faceZones().findZoneID(setName) >= 0)
    {
        return topoSetSource::sourceType::FACEZONE_SOURCE;
    }

    if (mesh().pointZones().findZoneID(setName) >= 0)
    {
        return topoSetSource::sourceType::POINTZONE_SOURCE;
    }

    return topoSetSource::sourceType::UNKNOWN_SOURCE;
}


Foam::topoSetSource::sourceType
Foam::expressions::fvExprDriver::topoSourceType(const word& setName) const
{
    auto setType = topoZoneType(setName);

    if (topoSetSource::sourceType::UNKNOWN_SOURCE == setType)
    {
        setType = topoSetType(setName);
    }

    return setType;
}



bool Foam::expressions::fvExprDriver::isCellSet(const word& setName) const
{
    return
    (
        topoSetSource::sourceType::CELLSET_SOURCE
     == topoSetType(setName)
    );
}


bool Foam::expressions::fvExprDriver::isFaceSet(const word& setName) const
{
    return
    (
        topoSetSource::sourceType::FACESET_SOURCE
     == topoSetType(setName)
    );
}


bool Foam::expressions::fvExprDriver::isPointSet(const word& setName) const
{
    return
    (
        topoSetSource::sourceType::POINTSET_SOURCE
     == topoSetType(setName)
    );
}


bool Foam::expressions::fvExprDriver::isCellZone(const word& name) const
{
    return (mesh().cellZones().findZoneID(name) >= 0);
}


bool Foam::expressions::fvExprDriver::isFaceZone(const word& name) const
{
    return (mesh().faceZones().findZoneID(name) >= 0);
}


bool Foam::expressions::fvExprDriver::isPointZone(const word& name) const
{
    return (mesh().pointZones().findZoneID(name) >= 0);
}


const Foam::expressions::exprResult&
Foam::expressions::fvExprDriver::lookupGlobal
(
    const word& name
) const
{
    return exprResultGlobals::New(this->mesh()).get(name, globalScopes_);
}


bool Foam::expressions::fvExprDriver::hasDataToWrite() const
{
    return (!storedVariables_.empty() || !delayedVariables_.empty());
}


void Foam::expressions::fvExprDriver::getData
(
    const dictionary& dict
)
{
    dict.readIfPresent("storedVariables", storedVariables_);
}


void Foam::expressions::fvExprDriver::prepareData
(
    dictionary& dict
) const
{
    auto& driver = const_cast<fvExprDriver&>(*this);

    (void)driver.update();

    if (storedVariables_.size())
    {
        driver.updateSpecialVariables(true);

        dict.add("storedVariables", storedVariables_);
    }
}


// ************************************************************************* //

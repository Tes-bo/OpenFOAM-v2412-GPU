/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "writeObjects.H"
#include "Time.H"
#include "polyMesh.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeObjects, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeObjects,
        dictionary
    );
}
}

const Foam::Enum
<
    Foam::functionObjects::writeObjects::writeOption
>
Foam::functionObjects::writeObjects::writeOptionNames_
({
    { writeOption::NO_WRITE, "noWrite" },
    { writeOption::AUTO_WRITE, "autoWrite" },
    { writeOption::ANY_WRITE, "anyWrite" },
    { writeOption::LOG, "log" },
});

const Foam::objectRegistry& setRegistry
(
    const Foam::Time& runTime,
    const Foam::dictionary& dict
)
{
    const Foam::word regionName =
        dict.getOrDefault("region", Foam::polyMesh::defaultRegion);

    if (regionName == "__TIME__")
    {
        return runTime;
    }

    return runTime.lookupObject<Foam::objectRegistry>(regionName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeObjects::writeObjects
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    obr_(setRegistry(runTime, dict)),
    writeOption_(ANY_WRITE),
    objectNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeObjects::read(const dictionary& dict)
{
    if (!functionObject::read(dict))
    {
        return false;
    }

    if (dict.found("field"))
    {
        objectNames_.resize(1);
        dict.readEntry("field", objectNames_.first());
    }
    else if (dict.found("fields"))
    {
        dict.readEntry("fields", objectNames_);
    }
    else
    {
        dict.readEntry("objects", objectNames_);
    }

    writeOption_ = writeOptionNames_.getOrDefault
    (
        "writeOption",
        dict,
        writeOption::ANY_WRITE
    );

    return true;
}


bool Foam::functionObjects::writeObjects::execute()
{
    return true;
}


bool Foam::functionObjects::writeObjects::write()
{
    if (writeOption_ == writeOption::LOG)
    {
        const auto& classes = obr_.classes();

        Log << "Registered objects:\n";

        forAllConstIters(classes, classInfo)
        {
            const word& className = classInfo.key();
            const wordHashSet& objectSet = classInfo();

            Log << "    " << className << ":\n";

            for (const auto& objectName : objectSet)
            {
                Log << "        " << objectName << "\n";
            }
            Log << nl;
        }

        Log << endl;

        return true;
    }

    Log << type() << " " << name() << " write:" << nl;

    if (!obr_.time().writeTime())
    {
        obr_.time().writeTimeDict();
    }

    // Get selection
    const wordList selectedNames(obr_.sortedNames<regIOobject>(objectNames_));

    // Warning if anything was missed
    bitSet missed(objectNames_.size());

    label index = 0;
    for (const wordRe& select : objectNames_)
    {
        if (!ListOps::found(selectedNames, select))
        {
            missed.set(index);
        }
        ++index;
    }

    if (missed.any())
    {
        WarningInFunction
            << "No corresponding selection for "
            << flatOutput(subset(missed, objectNames_)) << nl
            << "Available objects in database:"
            << nl << obr_.sortedToc()
            << endl;
    }

    for (const word& objName : selectedNames)
    {
        regIOobject& obj = obr_.lookupObjectRef<regIOobject>(objName);

        switch (writeOption_)
        {
            case writeOption::NO_WRITE:
            {
                if (obj.writeOpt() != IOobject::NO_WRITE)
                {
                    continue;
                }

                break;
            }
            case writeOption::AUTO_WRITE:
            {
                if (obj.writeOpt() != IOobject::AUTO_WRITE)
                {
                    continue;
                }

                break;
            }
            case writeOption::ANY_WRITE:
            {
                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown writeOption "
                    << writeOptionNames_[writeOption_]
                    << ". Valid writeOption types are "
                    << writeOptionNames_
                    << exit(FatalError);

                continue;
                break;
            }
        }

        if
        (
            obj.writeOpt() == IOobject::AUTO_WRITE
         && obr_.time().writeTime()
        )
        {
            Log << "    automatically written object " << obj.name() << endl;
        }
        else
        {
            // TBD:
            // If the object is a temporary field expression wrap with tmp<...>

            // if (obj.db().is_cacheTemporaryObject(obj))
            // {
            //     const word oldName(obj.name());
            //     obj.IOobject::rename("tmp<" + oldName + ">");
            //
            //     Log << "    writing object " << obj.name() << endl;
            //     obj.write();
            //     obj.IOobject::rename(oldName);
            // }
            // else
            {
                Log << "    writing object " << obj.name() << endl;
                obj.write();
            }
        }
    }

    return true;
}


// ************************************************************************* //

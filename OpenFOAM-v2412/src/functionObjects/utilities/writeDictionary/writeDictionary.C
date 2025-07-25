/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "writeDictionary.H"
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeDictionary, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeDictionary,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeDictionary::writeHeader()
{
    if (firstChange_)
    {
        firstChange_ = false;

        Info<< type() << ' ' << name() << " write:" << nl << nl;
        IOobject::writeDivider(Info) << endl;
    }
}


void Foam::functionObjects::writeDictionary::checkDictionary
(
    const dictionary& dict,
    const label dicti
)
{
    const auto digest(dict.digest());

    if (digests_[dicti] != digest)
    {
        digests_[dicti] = digest;

        writeHeader();
        dict.writeEntry(Info);
        Info<< nl;
        IOobject::writeDivider(Info) << endl;
    }
}


bool Foam::functionObjects::writeDictionary::tryDirectory
(
    const word& location,
    const label dicti
)
{
    IOobject dictIO
    (
        dictNames_[dicti],
        location,
        obr_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (dictIO.typeHeaderOk<IOdictionary>(true))
    {
        checkDictionary(IOdictionary(dictIO), dicti);

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeDictionary::writeDictionary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    dictNames_(),
    digests_(),
    firstChange_(true)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeDictionary::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    // Make unique
    dictNames_ = wordHashSet(dict.get<wordList>("dictNames")).sortedToc();

    digests_.resize(dictNames_.size());
    digests_ = SHA1Digest();

    Info<< type() << ' ' << name() << ": monitoring dictionaries:" << nl;
    for (const word& dictName : dictNames_)
    {
        Info<< "    " << dictName << nl;
    }
    if (dictNames_.empty())
    {
        Info<< "    none" << nl;
    }
    Info<< endl;

    return true;
}


bool Foam::functionObjects::writeDictionary::performCheck()
{
    // Restart reporting cycle
    firstChange_ = true;

    forAll(dictNames_, dicti)
    {
        // Also search parent (eg, Time) as too
        const auto* dictptr =
            obr_.cfindObject<IOdictionary>(dictNames_[dicti], true);

        if (dictptr)
        {
            checkDictionary(*dictptr, dicti);
        }
        else if (dictNames_[dicti] == Time::controlDictName)
        {
            // Slight hack. controlDict an unwatchedIOdictionary
            // (not IOdictionary) and not registered on Time either
            // - grab directly from Time
            checkDictionary(obr_.time().controlDict(), dicti);
        }
        else
        {
            const bool ok
            (
                tryDirectory(obr_.time().timeName(), dicti)
             || tryDirectory(obr_.time().constant(), dicti)
             || tryDirectory(obr_.time().system(), dicti)
            );

            if (!ok)
            {
                writeHeader();

                Info<< "    Unable to locate dictionary "
                    << dictNames_[dicti] << nl << nl;

                IOobject::writeDivider(Info) << endl;
            }
        }
    }

    return true;
}


bool Foam::functionObjects::writeDictionary::execute()
{
    return true;
}


bool Foam::functionObjects::writeDictionary::write()
{
    performCheck();
    return true;
}


// ************************************************************************* //

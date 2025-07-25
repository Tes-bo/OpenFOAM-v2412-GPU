/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

Description
    Basic tests of IOobjectList

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "volFields.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "hashedWordList.H"
#include "labelIOList.H"
#include "scalarIOList.H"

using namespace Foam;


void report(const UPtrList<const IOobject>& objects)
{
    Info<< "name/type:" << nl
        << objects.size() << nl << '(' << nl;

    for (const IOobject& io : objects)
    {
        Info<< "  " << io.name() << " : " << io.headerClassName() << nl;
    }

    Info<< ')' << nl << endl;
}


template<class Type>
void report(const UPtrList<const IOobject>& objects)
{
    Info<< "name/type:" << nl
        << objects.size() << nl << '(' << nl;

    for (const IOobject& io : objects)
    {
        Info<< "  " << io.name() << " : " << io.headerClassName()
            << (io.isHeaderClass<Type>() ? " is " : " is not ")
            << Type::typeName << nl;
    }

    Info<< ')' << nl << endl;
}


void report(const IOobjectList& objects)
{
    Info<< "Names: " << flatOutput(objects.sortedNames()) << nl
        << "Objects: " << objects << nl
        << "----" << nl;
}


void reportDetail(const IOobjectList& objects)
{
    Info<<"Details:" << nl;

    for (const word& key : objects.sortedNames())
    {
        // Canonical method name (NOV-2018)
        const IOobject* io = objects.findObject(key);

        label count = 0;

        // Test deprecated alternatives
        {
            // (const char*)
            IOobject* ptr = objects.lookup("SomeNonExistentName");
            if (ptr) ++count;
        }
        {
            // (const word&)
            IOobject* ptr = objects.lookup(key);
            if (ptr) ++count;
        }

        Info<< key << " (" << io->headerClassName()
            << ") = addr " << name(io) << nl;

        if (count != 1)
        {
            Warning
                << key << " had incorrect lookup?" << nl;
        }
    }

    Info<<"====" << nl;
}


void printFound(const IOobject* ptr)
{
    Info<< (ptr ? "found" : "not found") << nl;
}


void findObjectTest(const IOobjectList& objs)
{
    Info<< "Test findObject()" << nl << nl;

    const int oldDebug = IOobject::debug;
    IOobject::debug = 1;

    {
        Info<< "cfindObject(U)" << nl;
        const IOobject* io = objs.cfindObject("U");
        printFound(io);
    }

    {
        Info<< "getObject(U)" << nl;
        IOobject* io = objs.getObject("U");
        printFound(io);
    }

    {
        Info<< "cfindObject<void>(U)" << nl;
        const IOobject* io = objs.cfindObject<void>("U");
        printFound(io);
    }

    {
        Info<< "cfindObject<volScalarField>(U)" << nl;
        const IOobject* io = objs.cfindObject<volScalarField>("U");
        printFound(io);
    }

    {
        Info<< "cfindObject<volVectorField>(U)" << nl;
        const IOobject* io = objs.cfindObject<volVectorField>("U");
        printFound(io);
    }

    Info<< nl;

    IOobject::debug = oldDebug;
}


template<class Type>
void filterTest(const IOobjectList& objs, const wordRe& re)
{
    Info<< "Filter = " << re << nl;

    const word& typeName = Type::typeName;

    Info<< "    <" << typeName <<">(" << re << ") : "
        << objs.count<Type>(re) << nl
        << "    (" << typeName << "::typeName, " << re << ") : "
        << objs.count(typeName, re) << nl;

    Info<< "    <" << typeName << ">(" << re << ") : "
        << flatOutput(objs.sortedNames<Type>(re)) << nl
        // << flatOutput(objs.names<Type>(re)) << nl
        << "    (" << typeName << "::typeName, " << re << ") : "
        << flatOutput(objs.sortedNames(typeName, re)) << nl
        //<< flatOutput(objs.names(typeName, re)) << nl
        ;


    wordRe reClass("vol.*Field", wordRe::REGEX);
    wordRe re2(re, wordRe::REGEX_ICASE);

    Info<< "General" << nl
        << "    <void>(" << re << ") : "
        << flatOutput(objs.sortedNames<void>(re)) << nl
        << "    (" << reClass << ", " << re2 <<" ignore-case) : "
        << flatOutput(objs.sortedNames(reClass, re2)) << nl
        ;

    Info<< nl;
}


void registryTests(const IOobjectList& objs)
{
    Info<< nl << "IOobjectList " << flatOutput(objs.sortedNames()) << nl;

    Info<< "count" << nl
        << "    <void>()    : " << objs.count<void>() << nl
        << "    <labelList>()   : " << objs.count<labelIOList>() << nl
        << "    <scalarList>()   : " << objs.count<scalarIOList>() << nl
        << nl;
    Info<< "    <volScalarField>()    : "
        << objs.count<volScalarField>() << nl
        << "    (volScalarField::typeName) : "
        << objs.count(volScalarField::typeName) << nl;
    Info<< "    <volVectorField>()    : "
        << objs.count<volVectorField>() << nl
        << "    (volVectorField::typeName) : "
        << objs.count(volVectorField::typeName) << nl;


    Info<< nl << "Filter on names:" << nl;
    filterTest<volScalarField>(objs, wordRe("[p-z].*", wordRe::DETECT));

    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // argList::noParallel();
    argList::addOption
    (
        "filter",
        "wordRes",
        "filter keys with names or regexs"
    );
    argList::addBoolOption
    (
        "merge",
        "test merging lists (requires -filter)"
    );

    // timeSelector::addOptions();
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"
    #include "createTime.H"

    wordRes matcher;
    if (args.readListIfPresent<wordRe>("filter", matcher))
    {
        Info<<"limit names: " << matcher << nl;
    }

    if (args.found("merge") && matcher.empty())
    {
        FatalError
            << nl << "The -merge test also requires -filter" << nl
            << exit(FatalError);
    }


    const hashedWordList subsetTypes
    {
        volScalarField::typeName,
        volScalarField::Internal::typeName,
        volVectorField::typeName,
    };


    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        // Objects at this time
        IOobjectList objects(runTime, runTime.timeName());
        HashTable<wordHashSet> classes =
        (
            matcher.size()
          ? objects.classes(matcher)
          : objects.classes()
        );

        Info<< "Time: " << runTime.timeName() << nl;

        report(objects);
        report(objects.csorted());

        report(objects.csorted<volScalarField>());
        report(objects.csorted<volVectorField>());

        // Extra checks
        report<volScalarField>(objects.csorted<volScalarField>());
        report<volScalarField>(objects.csorted<volVectorField>());


        findObjectTest(objects);

        classes.filterKeys(subsetTypes);
        Info<<"only retain: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;

        classes = objects.classes();
        classes.erase(subsetTypes);
        Info<<"remove: " << flatOutput(subsetTypes) << nl;
        Info<<"Pruned: " << classes << nl;

        registryTests(objects);

        // On last time
        if (timeI == timeDirs.size()-1)
        {
            if (args.found("merge"))
            {
                Info<< nl << "Test merge" << nl;
            }
            else
            {
                continue;
            }

            IOobjectList other(runTime, runTime.timeName());

            Info<< "==original==" << nl; reportDetail(objects);

            objects.filterKeys(matcher);

            Info<< "==target==" << nl; reportDetail(objects);
            Info<< "==source==" << nl; reportDetail(other);

            objects.merge(other);
            Info<< nl << "After merge" << nl;

            Info<< "==target==" << nl; reportDetail(objects);
            Info<< "==source==" << nl; reportDetail(other);

            Info<< nl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

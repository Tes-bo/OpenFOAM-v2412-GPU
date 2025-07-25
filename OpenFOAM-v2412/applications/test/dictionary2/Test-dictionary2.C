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

Application
    Test-dictionary2

Description
    Test dictionary insertion and some reading functionality.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"
#include "dictionary.H"
#include "SpanStream.H"
#include "ops.H"
#include "scalarRange.H"
#include "stringOps.H"

using namespace Foam;

void entryInfo(entry* e)
{
    if (e)
    {
        Info<<"added "
            << e->keyword() << ": " << typeid(e).name() << nl;
    }
}


// Try with readScalar
scalar try_readScalar(const dictionary& dict, const word& k)
{
    scalar val(-GREAT);

    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOerr = FatalIOError.throwing(true);

    try
    {
        val = readScalar(dict.lookup(k));
        Info<< "readScalar(" << k << ") = " << val << nl;
    }
    catch (const Foam::IOerror& err)
    {
        Info<< "readScalar(" << k << ") Caught FatalIOError "
            << err << nl << endl;
    }
    catch (const Foam::error& err)
    {
        Info<< "readScalar(" << k << ") Caught FatalError "
            << err << nl << endl;
    }
    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOerr);

    return val;
}


// Try with get<scalar>
scalar try_getScalar(const dictionary& dict, const word& k)
{
    scalar val(-GREAT);

    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOerr = FatalIOError.throwing(true);

    try
    {
        val = dict.get<scalar>(k);
        Info<< "get<scalar>(" << k << ") = " << val << nl;
    }
    catch (const Foam::IOerror& err)
    {
        Info<< "get<scalar>(" << k << ") Caught FatalIOError "
            << err << nl << endl;
    }
    catch (const Foam::error& err)
    {
        Info<< "get<scalar>(" << k << ") Caught FatalError "
            << err << nl << endl;
    }
    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOerr);

    return val;
}


// Try with getCheck<scalar>
template<class Predicate>
scalar try_getCheckScalar
(
    const dictionary& dict,
    const word& k,
    const Predicate& pred
)
{
    scalar val(-GREAT);

    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOerr = FatalIOError.throwing(true);

    try
    {
        val = dict.getCheck<scalar>(k, pred);
        Info<< "getCheck<scalar>(" << k << ") = " << val << nl;
    }
    catch (const Foam::IOerror& err)
    {
        Info<< "getCheck<scalar>(" << k << ") Caught FatalIOError "
            << err << nl << endl;
    }
    catch (const Foam::error& err)
    {
        Info<< "getCheck<scalar>(" << k << ") Caught FatalError "
            << err << nl << endl;
    }
    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOerr);

    return val;
}


// Try with *entry (from findEntry) and get<scalar>
scalar try_getScalar(const entry* eptr, const word& k)
{
    scalar val(-GREAT);

    if (!eptr)
    {
        Info<< "No entry" << k << nl;
        return val;
    }

    const bool oldThrowingError = FatalError.throwing(true);
    const bool oldThrowingIOerr = FatalIOError.throwing(true);

    try
    {
        val = eptr->get<scalar>();
        Info<< "entry get<scalar>(" << k << ") = " << val << nl;
    }
    catch (const Foam::IOerror& err)
    {
        Info<< "entry get<scalar>(" << k << ") Caught FatalIOError "
            << err << nl << endl;
    }
    catch (const Foam::error& err)
    {
        Info<< "entry get<scalar>(" << k << ") Caught FatalError "
            << err << nl << endl;
    }
    FatalError.throwing(oldThrowingError);
    FatalIOError.throwing(oldThrowingIOerr);

    return val;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();

    argList args(argc, argv);

    dictionary dict1;
    for (label i=0; i<5; ++i)
    {
        dictionary tmpdict;

        {
            entry* e = dict1.add
            (
                word::printf("entry%d", i),
                string("entry" + Foam::name(i))
            );
            entryInfo(e);
        }

        {
            entry* e = tmpdict.add
            (
                word::printf("subentry%d", i),
                string("subentry" + Foam::name(i))
            );
            entryInfo(e);
        }

        {
            entry* e = dict1.add
            (
                word::printf("dict%d", i),
                tmpdict
            );
            entryInfo(e);
        }
    }

    // Insert new dictionary or merge into existing one
    for (auto k : { "dict1", "dict10" })
    {
        const word key(k);
        entry* e = dict1.add(key, dictionary(), true);

        if (e && e->isDict())
        {
            e->dict().add(word("sub1" + key), 10);
            e->dict().add(word("sub2" + key), 20);
            e->dict().add(word("sub3" + key), 30);
            e->dict().add(word("sub4" + key), 40);
            e->dict().add(word("sub5" + key), 50);
            e->dict().add(word("sub6" + key), 60);
        }
    }


    // overwrite existing
    {
        dict1.set("entry3", 1000); // overwrite
        entry* e = dict1.set(word("dict3"), 1000); // overwrite
        entryInfo(e);
    }

    // merge into existing dictionary: returns pointer to existing dict
    {
        dictionary tmpdict;
        tmpdict.add(word("something"), 3.14159);

        entry* e = dict1.add(word("dict4"), tmpdict, true);  // merge
        entryInfo(e);

        if (e) Info<< nl << "=> " << *e << nl;

        tmpdict.clear();
        tmpdict.add(word("other"), 2.718281);

        dict1.add(word("dict1"), tmpdict, true);  // merge
    }

    Info<< nl << "dictionary" << nl << nl;
    dict1.write(Info, false);


    {
        dict1.foundCompat
        (
            "newEntry", {{"entry1", 1612}, {"entry15", 1606}}
        );
        dict1.foundCompat
        (
            "newEntry", {{"entry15", 1612}, {"entry2", 1606}}
        );
        dict1.foundCompat
        (
            "newEntry", {{"entry3", 240}, {"entry2", 1606}}
        );

        // And some success
        dict1.foundCompat
        (
            "entry4", {{"none", 240}, {"entry2", 1606}}
        );
    }

    {
        label lval = readLabel
        (
            dict1.lookupCompat
            (
                "entry400", {{"none", 240}, {"entry3", 1606}}
            )
        );
        Info<< "int value: " << lval << nl;
    }


    // Could have different dictionary names and different entries names etc.
    // Quite ugly!
    {
        scalar sval = readScalar
        (
            dict1.csearchCompat
            (
                "newdictName", {{"dict4", 1706}, {"dict1", 1606}}
            )
            .dict()
            .lookupCompat
            (
                "newval", {{"something", 1606}, {"other", 1612}}
            )
        );

        Info<< "scalar value: " << sval << nl;

        sval = readScalar
        (
            dict1.csearchCompat
            (
                "newdictName", {{"dict1", 1606}, {"dict4", 1706}}
            )
            .dict()
            .lookupCompat
            (
                "newval", {{"something", 1606}, {"other", 1612}}
            )
        );

        Info<< "scalar value = " << sval << nl;
    }


    Info<< nl << "dictionary" << nl << nl;
    dict1.write(Info, false);


    {
        Info<< nl << "Test reading good/bad/empty scalar entries" << nl;
        dictionary dict2
        (
            ICharStream
            (
                "good 3.14159;\n"
                "negative -3.14159;\n"
                "neg2 -3.14159;\n"
                "empty;\n"
                // "bad  text;\n"            // always fails
                // "bad  3.14159 1234;\n"    // fails for readScalar
            )()
        );
        dict2.write(Info);


        // With readScalar
        {
            Info<< nl << "Test some bad input with readScalar()" << nl;

            try_readScalar(dict2, "good");
            // try_readScalar(dict2, "bad");
            try_readScalar(dict2, "empty");
        }


        // With get<scalar>
        {
            Info<< nl << "Test some bad input with get<scalar>()" << nl;

            try_getScalar(dict2, "good");
            // try_getScalar(dict2, "bad");
            try_getScalar(dict2, "empty");
        }


        // With getCheck<scalar>
        {
            Info<< nl << "Test some input with getCheck<scalar>()" << nl;

            try_getCheckScalar(dict2, "good", scalarRange::gt0());
            try_getCheckScalar(dict2, "negative", scalarRange::gt0());
            try_getCheckScalar(dict2, "neg2", scalarRange::gt0());

            try_getCheckScalar(dict2, "good", greaterOp1<scalar>(0));
            try_getCheckScalar(dict2, "negative", greaterOp1<scalar>(0));

            Info<< nl << "with lambda" << nl;
            try_getCheckScalar
            (
                dict2,
                "good",
                [](const scalar x) { return x > 0; }
            );
        }

        // With findEntry and get<scalar>
        {
            Info<< nl
                << "Test some bad input with findEntry + get<scalar>()" << nl;

            try_getScalar(dict2.findEntry("good"), "good");
            // try_getScalar(dict2.findEntry("bad"), "bad");
            try_getScalar(dict2.findEntry("empty"), "empty");
        }

        #ifdef COMPAT_OPENFOAM_ORG
        {
            Info<< nl
                << "Test openfoam.org compatibility" << nl;

            dictionary mydict
            (
                ICharStream
                (
                    "scalar 3.14159;\n"
                    "label 10;\n"
                )()
            );

            Info<<"get<scalar> : " << mydict.get<scalar>("scalar") << nl;
            Info<<"get<label> : " << mydict.get<label>("label") << nl;

            Info<<"lookup<scalar> : " << mydict.lookup<scalar>("scalar") << nl;
            Info<<"lookup<label> : " << mydict.lookup<label>("label") << nl;
        }
        #else
        Info<< "No openfoam.org compatibility methods" << nl;
        #endif
    }


    Info<< "\nDone\n" << endl;

    return 0;
}


// ************************************************************************* //

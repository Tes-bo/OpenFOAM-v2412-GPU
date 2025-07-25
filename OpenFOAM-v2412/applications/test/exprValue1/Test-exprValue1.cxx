/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Application
    Test-exprValue1

Description
    Test low-level polymorphic value container (exprValue)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "ITstream.H"
#include "OTstream.H"
#include "SpanStream.H"
#include "exprValue.H"
#include "Pstream.H"

using namespace Foam;

void printInfo(const expressions::exprValue& val)
{
    Pout<< "Boxed type:" << int(val.typeCode())
        << " (" << val.valueTypeName() << ") good:"
        << val.good() << " => " << val << nl;
}


void write_read(const expressions::exprValue& val)
{
    OCharStream os;
    os << val;

    ISpanStream is(os.view());
    expressions::exprValue val2;
    is >> val2;

    Pout<< "wrote " << os.count() << " chars: " << os.str() << nl;

    printInfo(val);
    printInfo(val2);
    Pout<< "====" << nl;
}


tokenList tokens_of(const expressions::exprValue& val)
{
    OTstream toks;
    toks << val;

    Pout<< "val with tokens: ";
    toks.writeList(Pout, 0) << nl;

    for (const auto& t : toks)
    {
        Pout<< "  " << t.info() << nl;
    }
    Pout<< nl;

    return toks;
}


expressions::exprValue tryParse(const std::string& str)
{
    expressions::exprValue val, val2;

    ITstream is(str);

    const bool ok = val.read(is);

    Info<< "read " << Foam::name(val.typeCode()) << " from " << str;

    if (ok)
    {
        Info<< " trailing tokens:" << is.nRemainingTokens() << nl
            << "value: " << val << nl;
    }
    else
    {
        Info<< " FAILED" << nl;
    }

    if (ok)
    {
        Info<< "Direct from string: ";
        if (expressions::exprValue::read(str, val2))
        {
            Info<< "OK" << nl;
        }
        else
        {
            Info<< "NOK" << nl;
        }
    }
    return val;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    // Aborts
    // expressions::exprValue value(std::string(""));

    // Regular broadcast doesn't work
    Info<< "exprValue"
        << " sizeof:" << sizeof(expressions::exprValue)
        << " contiguous:" << is_contiguous<expressions::exprValue>::value
        << nl << nl;

    {
        expressions::exprValue value;
        tokenList toks;

        Info<< "exprValue"
            << " sizeof:" << value.size_bytes()
            << nl << nl;

        // Info<< "value: " << value << nl;

        // Nothing
        printInfo(value);
        toks = tokens_of(value);
        write_read(value);

        value.set(scalar(100));
        printInfo(value); write_read(value); toks = tokens_of(value);

        value.set(scalar(100.01));
        printInfo(value); write_read(value); toks = tokens_of(value);

        value.set(vector(1,2,3));
        printInfo(value); write_read(value); toks = tokens_of(value);

        value = vector(4,5,6);
        printInfo(value); write_read(value); toks = tokens_of(value);

        value = Zero;
        printInfo(value); write_read(value); toks = tokens_of(value);

        value.clear();
        printInfo(value); write_read(value); toks = tokens_of(value);

        value.set<bool>(true);

        printInfo(value);
        printInfo(value); write_read(value); toks = tokens_of(value);

        if (UPstream::parRun())
        {
            Info<< "Before broadcast" << nl;
        }

        if (UPstream::master())
        {
            value = 100 * vector(1,0,0);
        }
        printInfo(value);

        expressions::exprValue oldValue(value);


        if (true)
        {
            // Broadcast with serialization
            Pstream::broadcast(value);
        }
        else
        {
            // Broadcast manually
            UPstream::broadcast
            (
                value.data_bytes(),
                value.size_bytes(),
                UPstream::worldComm
            );
        }

        Pout<< "same values: " << (oldValue == value) << nl;

        // Can put them in a Map and write, but reading will not work...
        Map<expressions::exprValue> map;
        map(1) = value;
        Info<< "map: " << map << nl;

        if (UPstream::parRun())
        {
            Info<< "After broadcast" << nl;
            printInfo(value);
        }
    }


    {
        Info<< nl << "Test parsing" << nl << nl;

        for
        (
            const auto& input :
            stringList
            ({
                "()",  // bad
                "( 1 2  ",  // also bad
                "(  ",  // really bad
                "(1 16 12)",
                "(1 bad)",
                "(5)",
                "1.2345",
                "5.678 trailing",
                "true",
                "false",
                " 1 ",
                "  yes no "
            })
        )
        {
            (void) tryParse(input);
        }
    }

    return 0;
}


// ************************************************************************* //

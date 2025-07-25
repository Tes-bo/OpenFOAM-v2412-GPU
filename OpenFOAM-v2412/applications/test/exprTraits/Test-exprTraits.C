/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Basic tests of expression traits

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "ITstream.H"
#include "uLabel.H"

#include "exprTraits.H"
#include "error.H"
#include "stringList.H"
#include "exprScanToken.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void printTraits()
{
    const auto typeCode = exprTypeTraits<Type>::value;

    Info<< "Type '" << pTraits<Type>::typeName
        << "' = code:" << int(typeCode)
        << " rank:" << exprTypeTraits<Type>::rank
        << " cmpt:" << exprTypeTraits<Type>::nComponents
        << " name:" << exprTypeTraits<Type>::name;

    if (pTraits<Type>::typeName != word(exprTypeTraits<Type>::name))
    {
        Info<< " (UNSUPPORTED)";
    }

    Info << endl;
}


void print(const expressions::scanToken& tok)
{
    Info<< "    type:" << int(tok.type_);
    if (tok.is_pointer())
    {
        Info<< " ptr:" << Foam::name(tok.name_);
    }
    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    Info<< nl << "Traits:" << nl;

    printTraits<word>();
    printTraits<string>();
    printTraits<bool>();
    printTraits<label>();
    printTraits<scalar>();
    printTraits<complex>();
    printTraits<vector>();
    printTraits<tensor>();
    printTraits<symmTensor>();
    printTraits<sphericalTensor>();

    const auto getName = nameOp<expressions::valueTypeCode>();

    Info<< nl;

    Info<< "Name of typeCode: "
        << Foam::name(expressions::valueTypeCode::type_scalar) << nl;

    Info<< "Name of typeCode: "
        << getName(expressions::valueTypeCode::type_bool) << nl;

    {
        expressions::scanToken tok(expressions::scanToken::null());
        expressions::scanToken tok2(expressions::scanToken::null());

        Info<< nl << "sizeof(scanToken): "
            << sizeof(tok) << nl;

        print(tok);
        print(tok2);

        tok.setWord("hello");

        print(tok);

        tok2 = tok;
        print(tok2);

        tok2.destroy();

        print(tok);  // Not a leak, but old rubbish
        print(tok2);
    }

    Info<< nl << "Done" << nl;
    return 0;
}


// ************************************************************************* //

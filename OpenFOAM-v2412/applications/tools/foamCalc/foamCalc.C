/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Application
    foamCalc

Group
    grpTools

Description
    A simple expression calculator using OpenFOAM string evaluation.
    Multiple arguments will be concatenated together.

Usage
    \b foamCalc [OPTION]

    Options:
      - \par -precision int
        Output with specified precision

      - \par -rules
        Print parser rules and exit

      - \par -tokens
        Print token names  and exit

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "exprString.H"
#include "stringOps.H"
#include "fieldExprParser.H"

using namespace Foam;

template<class Parser>
void printInformation
(
    Ostream& os,
    const word& name,
    const bool printNames,
    const bool printRules
)
{
    if (printNames)
    {
        os << nl << "# Tokens for " << name << nl;
        Parser::printTokenNames(os);
    }

    if (printRules)
    {
        os << nl << "# Rules for " << name << nl;
        Parser::printRules(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::setAdvanced("case");  // Hide -case : has no meaning here

    argList::addNote
    (
        "A simple expression calculator using OpenFOAM string evaluation.\n"
        "Multiple arguments will be concatenated together."
    );

    argList::addBoolOption("rules", "Print parser rules and exit");
    argList::addBoolOption("tokens", "Print token names and exit");
    argList::addOption("precision", "int", "Output with specified precision");
    argList::addOption("size", "int", "Field output width (default: 1)");

    // Flag arguments as optional so that -rules and -tokens works
    argList::noMandatoryArgs();
    argList::addArgument
    (
        "expression",
        "The expression to evaluate"
    );

    argList args(argc, argv);

    const bool printRules = args.found("rules");
    const bool printNames = args.found("tokens");
    const label fieldWidth = args.getOrDefault<label>("size", 1);

    if (printNames || printRules)
    {
        printInformation<expressions::fieldExpr::parser>
        (
            Info,
            "field",
            printNames,
            printRules
        );

        return 0;
    }

    {
        const int prec = args.getOrDefault<int>("precision", 0u);
        if (prec > 0)
        {
            IOstream::defaultPrecision(prec);
            Sout.precision(prec);
        }
    }

    std::string expr;

    for (int argi=1; argi < args.size(); ++argi)
    {
        if (argi > 1) expr += ' ';
        expr += args[argi];
    }

    // Don't bother stripping C/C++, but do allow empty variables
    stringOps::inplaceExpand(expr, true);
    stringOps::inplaceTrim(expr);

    if (expr.empty())
    {
        InfoErr
            << "Error: no expression specified." << nl
            << "See '" << args.executable() << " -help' for usage" << nl
            << nl;

        return 1;
    }

    Info<< stringOps::evaluate(fieldWidth, expr).c_str() << nl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    foamExprParserInfo

Group
    grpTools

Description
    Display token names or rules for specified expression parsers.

    In the Lemon grammar, terminals (uppercase) are listed first.
    Non-terminals (lowercase) are listed second.

    The current OpenFOAM grammar short naming conventions:
    - svalue : scalar value
    - sfield : scalar field
    - lfield : logic field
    - vfield : vector field
    - tfield : tensor field
    - hfield : sphericalTensor field
    - yfield : symmTensor field
    .
    Prefixes: 's' (surface) or 'p' (point).
    For example, psfield for a point scalar field

Usage
    \b foamExprParserInfo [OPTION]

    Options:
      - \par -rules
        Print parser rules

      - \par -tokens
        Print token names (default)

      - \par -all
        Display information for all parsers

      - \par -field
        Field expression parser information

      - \par -patch
        Patch expression parser information

      - \par -volume
        Volume expression parser information

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOstreams.H"
#include "List.H"
#include "fieldExprParser.H"
#include "patchExprParser.H"
#include "volumeExprParser.H"

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
        "Display token names or rules for specified expression parsers.\n"
        "In the Lemon grammar, terminals (uppercase) are listed first.\n"
        "Non-terminals (lowercase) are listed second.\n \n"
        "The current OpenFOAM grammar short naming conventions:\n"
        "  * svalue : scalar value\n"
        "  * sfield : scalar field\n"
        "  * lfield : logic field\n"
        "  * vfield : vector field\n"
        "  * tfield : tensor field\n"
        "  * hfield : sphericalTensor field\n"
        "  * yfield : symmTensor field\n"
        " \n"
        "Prefixes: 's' (surface) or 'p' (point).\n"
        "Eg, psfield for a point scalar field\n"
    );

    argList::addBoolOption("rules", "Print parser rules");
    argList::addBoolOption("tokens", "Print token names (default)");

    argList::addBoolOption("all", "Display information for all parsers");
    argList::addBoolOption("field", "Field expression parser information");
    argList::addBoolOption("patch", "Patch expression parser information");
    argList::addBoolOption("volume", "Volume expression parser information");

    argList args(argc, argv);

    const bool all = args.found("all");
    const bool printRules = args.found("rules");
    const bool printNames = args.found("tokens") || !printRules;

    label count = 0;

    if (all || args.found("field"))
    {
        ++count;
        printInformation<expressions::fieldExpr::parser>
        (
            Info,
            "field",
            printNames,
            printRules
        );
    }

    if (all || args.found("patch"))
    {
        ++count;
        printInformation<expressions::patchExpr::parser>
        (
            Info,
            "patch",
            printNames,
            printRules
        );
    }

    if (all || args.found("volume"))
    {
        ++count;
        printInformation<expressions::volumeExpr::parser>
        (
            Info,
            "volume",
            printNames,
            printRules
        );
    }

    if (!count)
    {
        InfoErr
            << "Error: no parser selected." << nl
            << "See '" << args.executable() << " -help' for usage" << nl
            << nl;

        return 1;
    }

    return 0;
}


// ************************************************************************* //

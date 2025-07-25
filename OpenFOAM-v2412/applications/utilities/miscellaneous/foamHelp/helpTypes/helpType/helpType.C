/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2015 OpenFOAM Foundation
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

#include "helpType.H"
#include "doxygenXmlParser.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(helpType, 0);
    defineRunTimeSelectionTable(helpType, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::helpType::doxygenPath() const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    List<fileName> docDirs(docDict.lookup("doxyDocDirs"));

    label dirI = -1;
    forAll(docDirs, i)
    {
        if (isDir(docDirs[i].expand()))
        {
            dirI = i;
            break;
        }
    }

    if (dirI == -1)
    {
        Info<< "No Doxygen sources found under search paths: "
            << docDirs << endl;
        return fileName();
    }

    return docDirs[dirI];
}


void Foam::helpType::displayDocOptions
(
    const string& searchStr,
    const bool exactMatch,
    const word& ext
) const
{
    fileName doxyPath(doxygenPath());

    if (doxyPath.empty())
    {
        return;
    }

    Info<< "Found doxygen help in " << doxyPath.c_str() << nl << endl;

    doxygenXmlParser parser
    (
        doxyPath/"../DTAGS",
        "tagfile",
        searchStr,
        exactMatch,
        ext
    );

    if (debug)
    {
        Info<< parser;
    }

    Info<< "Valid types:" << nl << parser.sortedToc();
}


void Foam::helpType::displayDoc
(
    const word& className,
    const string& searchStr,
    const bool exactMatch,
    const word& ext
) const
{
    fileName doxyPath(doxygenPath());

    if (doxyPath.empty())
    {
        return;
    }

    Info<< "Found doxygen help in " << doxyPath.c_str() << nl << endl;

    string docBrowser = getEnv("FOAM_DOC_BROWSER");
    if (docBrowser.empty())
    {
        const dictionary& docDict =
            debug::controlDict().subDict("Documentation");
        docDict.readEntry("docBrowser", docBrowser);
    }

    doxygenXmlParser parser
    (
        doxyPath/"../DTAGS",
        "tagfile",
        searchStr,
        exactMatch,
        ext
    );

    if (debug)
    {
        Info<< parser;
    }

    if (parser.found(className))
    {
        fileName docFile
        (
            doxyPath/parser.subDict(className).get<fileName>("filename")
        );

        // can use FOAM_DOC_BROWSER='application file://%f' if required
        docBrowser.replaceAll("%f", docFile);

        fileName classDirectory
        (
            parser.subDict(className).get<fileName>("path")
        );
        const word classFile
        (
            parser.subDict(className).get<word>("name")
        );

        Info<< "Showing documentation for type " << className << nl << endl;

        Info<< "Source file: " << classDirectory.c_str() << classFile << nl
            << endl;

        Foam::system(docBrowser);
    }
    else
    {
        FatalErrorInFunction
            << "No help for type " << className << " found."
            << "  Valid options:" << parser.sortedToc()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpType::init()
{
    argList::addOption
    (
        "browse",
        "word",
        "Display documentation in browser"
    );
}


// ************************************************************************* //

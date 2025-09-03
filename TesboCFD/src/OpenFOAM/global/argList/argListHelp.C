/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "argList.H"
#include "foamVersion.H"
#include "stringOps.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Convert api number (YYMM) -> 20YY. Eg, 1906 -> 2019
static inline int apiYear()
{
    return 2000 + (foamVersion::api / 100);
}


// Footer for manpage
static inline void printManFooter()
{
    Info<< ".SH \"SEE ALSO\"" << nl
        << "Online documentation "
        << "https://www.openfoam.com/documentation/" << nl
        << ".SH COPYRIGHT" << nl
        << "Copyright \\(co 2018-" << apiYear() << " OpenCFD Ltd." << nl;
}


// Option output (manpage formatted)
static inline void printManOption
(
    const word& optName,
    const string& optArg,
    const string& optUsage
)
{
    Info<< ".TP\n\\fB\\-" << optName << "\\fR";

    if (optArg.size())
    {
        Info<< " \\fI" << optArg.c_str() << "\\fR";
    }
    Info<< nl;

    if (optUsage.size())
    {
        stringOps::writeWrapped(Info, optUsage, argList::usageMax, 0, true);
    }
    else
    {
        Info<< nl;
    }
}


// Bool option output (manpage formatted)
static inline void printManOption
(
    const word& optName,
    const string& optUsage
)
{
    printManOption(optName, string::null, optUsage);
}


// Option output (manpage formatted)
// - uses static HashTables to obtain values
static void printManOption(const word& optName)
{
    printManOption
    (
        optName,
        argList::validOptions.lookup(optName, string::null),
        argList::optionUsage.lookup(optName, string::null)
    );

    if (argList::validParOptions.contains(optName))
    {
        Info<< "\\fB[Parallel option]\\fR" << nl;
    }
}


// Wrapped output with initial start column
static void printOptionUsage
(
    std::string::size_type start,
    const string& str
)
{
    if (str.empty())
    {
        Info<< nl;
        return;
    }

    // Indent the first line. Min 2 spaces between option and usage
    if (start + 2 > argList::usageMin)
    {
        Info<< nl;
        start = 0;
    }
    while (start < argList::usageMin)
    {
        Info<< ' ';
        ++start;
    }

    stringOps::writeWrapped
    (
        Info,
        str,
        (argList::usageMax - argList::usageMin),
        argList::usageMin
    );
}


// Option output (usage formatted)
static inline void printOption
(
    const word& optName,
    const string& optArg,
    const string& optUsage
)
{
    Info<< "  -" << optName;

    // Length with leading '  -'
    std::string::size_type len = optName.size() + 3;

    if (optArg.size())
    {
        Info<< " <" << optArg.c_str() << '>';

        // Length with space between option/param and '<>'
        len += optArg.size() + 3;
    }

    printOptionUsage(len, optUsage);
}


// Bool option output (usage formatted)
static inline void printOption
(
    const word& optName,
    const string& optUsage
)
{
    printOption(optName, string::null, optUsage);
}


// Option output (usage formatted)
// - uses static HashTables to obtain values
static void printOption(const word& optName)
{
    printOption
    (
        optName,
        argList::validOptions.lookup(optName, string::null),
        argList::optionUsage.lookup(optName, string::null)
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::argList::printMan() const
{
    // .TH "<APPLICATION>" 1 "OpenFOAM-<version>" "source" "category"

    Info<< ".TH" << token::SPACE
        // All uppercase (returns a Foam::string) and thus also quoted
        << stringOps::upper(executable_) << token::SPACE
        << 1 << token::SPACE
        << token::DQUOTE << "OpenFOAM-v" << foamVersion::api << token::DQUOTE
        << token::SPACE
        << token::DQUOTE << "www.openfoam.com" << token::DQUOTE
        << token::SPACE
        << token::DQUOTE << "OpenFOAM Commands Manual" << token::DQUOTE
        << nl;


    // .SH NAME
    // <application> \- part of OpenFOAM (The Open Source CFD Toolbox).
    Info<< ".SH \"NAME\"" << nl
        << executable_
        << " \\- part of \\fBOpenFOAM\\fR (The Open Source CFD Toolbox)."
        << nl;


    // .SH SYNOPSIS
    // .B command [OPTIONS] ...

    Info<< ".SH \"SYNOPSIS\"" << nl
        << "\\fB" << executable_ << "\\fR [\\fIOPTIONS\\fR]";

    if (validArgs.size())
    {
        Info<< ' ';

        if (!argsMandatory_)
        {
            Info<< '[';
        }

        label i = 0;
        for (const std::string& argName : validArgs)
        {
            if (i++) Info<< ' ';
            Info << "\\fI" << argName.c_str() << "\\fR";
        }

        if (!argsMandatory_)
        {
            Info<< ']';
        }
    }
    Info<< nl;


    // .SH DESCRIPTION
    {
        Info<< ".SH \"DESCRIPTION\"" << nl;

        Info<< ".nf" << nl; // No fill lines

        if (notes.empty())
        {
            Info<< "No description available\n";
        }
        else
        {
            for (const std::string& note : notes)
            {
                if (note.empty())
                {
                    Info<< nl;
                }
                else
                {
                    stringOps::writeWrapped(Info, note, usageMax, 0, true);
                }
            }
        }
        Info<< ".fi" << nl; // Fill lines
    }


    // Arguments output
    if (validArgs.size())
    {
        // .SH "ARGUMENTS"
        Info<< ".SH \"ARGUMENTS\"" << nl;

        label argIndex = 0;
        for (const std::string& argName : validArgs)
        {
            ++argIndex;

            Info<< ".TP\n\\fI" << argName.c_str() << "\\fR";
            Info<< nl;

            const string& text =
                argList::argUsage.lookup(argIndex, string::null);

            if (text.size())
            {
                stringOps::writeWrapped(Info, text, usageMax, 0, true);
            }
            else
            {
                Info<< nl;
            }
        }
    }


    // .SH "OPTIONS"
    Info<< ".SH \"OPTIONS\"" << nl;

    for (const word& optName : validOptions.sortedToc())
    {
        // Normal (non-advanced) options
        if (!advancedOptions.contains(optName))
        {
            printManOption(optName);
        }
    }

    // Standard documentation/help options

    printManOption("doc", "Display documentation in browser");
    printManOption("help", "Display short help and exit");
    printManOption("help-full", "Display full help and exit");


    // .SS "ADVANCED OPTIONS"
    Info<< ".SS \"ADVANCED OPTIONS\"" << nl;

    for (const word& optName : validOptions.sortedToc())
    {
        // Advanced options
        if (advancedOptions.contains(optName))
        {
            printManOption(optName);
        }
    }

    printManOption("doc-source", "Display source code in browser");


    const bool hasCompat =
    (
        !argList::validOptionsCompat.empty()
     || !argList::ignoreOptionsCompat.empty()
    );

    if (hasCompat)
    {
        printManOption("help-compat", "Display compatibility options and exit");
    }

    printManOption("help-man", "Display full help (manpage format) and exit");
    printManOption("help-notes", "Display help notes (description) and exit");

    printManCompat();

    // Footer
    printManFooter();
}


void Foam::argList::printUsage(bool full) const
{
    Info<< "\nUsage: " << executable_ << " [OPTIONS]";

    if (validArgs.size())
    {
        Info<< ' ';

        if (!argsMandatory_)
        {
            Info<< '[';
        }

        label i = 0;
        for (const std::string& argName : validArgs)
        {
            if (i++) Info<< ' ';
            Info<< '<' << argName.c_str() << '>';
        }

        if (!argsMandatory_)
        {
            Info<< ']';
        }
    }
    Info<< nl;

    // Arguments output
    // Currently only if there is also usage information, but may wish to
    // change this to remind developers to add some description.
    if (validArgs.size() && argUsage.size())
    {
        Info<< "Arguments:\n";

        label argIndex = 0;
        for (const std::string& argName : validArgs)
        {
            ++argIndex;

            Info<< "  <" << argName.c_str() << '>';

            printOptionUsage
            (
                // Length with leading spaces and surround '<>'
                (argName.size() + 4),
                argList::argUsage.lookup(argIndex, string::null)
            );
        }
    }

    Info<< "Options:\n";

    for (const word& optName : validOptions.sortedToc())
    {
        // Suppress advanced options for regular -help.
        if (full || !advancedOptions.contains(optName))
        {
            printOption(optName);
        }
    }


    // Place documentation/help options at the end

    printOption("doc", "Display documentation in browser");
    if (full)
    {
        printOption("doc-source", "Display source code in browser");
    }
    printOption("help", "Display short help and exit");

    if
    (
        argList::validOptionsCompat.size()
      + argList::ignoreOptionsCompat.size()
    )
    {
        printOption("help-compat", "Display compatibility options and exit");
    }

    if (full)
    {
        printOption("help-man", "Display full help (manpage format) and exit");
        printOption("help-notes", "Display help notes (description) and exit");
    }

    printOption("help-full", "Display full help and exit");

    printNotes();

    Info<< nl;
    foamVersion::printBuildInfo(Info.stdStream(), true);
    Info<< endl;
}


void Foam::argList::printNotes() const
{
    // Output notes with automatic text wrapping
    if (!notes.empty())
    {
        Info<< nl;

        for (const std::string& note : notes)
        {
            if (note.empty())
            {
                Info<< nl;
            }
            else
            {
                stringOps::writeWrapped(Info, note, usageMax);
            }
        }
    }
}


void Foam::argList::printManCompat() const
{
    if
    (
        argList::validOptionsCompat.empty()
     && argList::ignoreOptionsCompat.empty()
    )
    {
        return;
    }


    // .SS "COMPATIBILITY OPTIONS"
    Info<< ".SS \"COMPATIBILITY OPTIONS\"" << nl;

    for (const word& k : argList::validOptionsCompat.sortedToc())
    {
        const auto& iter = *argList::validOptionsCompat.cfind(k);

        const word& optName = iter.first;
        const int until = abs(iter.second);

        Info<< ".TP\n\\fB\\-" << k
            << "\\fR (now \\fB\\-" << optName << "\\fR)" << nl;

        if (until)
        {
            Info<< "The option was last used in " << until << "." << nl;
        }
    }

    for (const word& k : argList::ignoreOptionsCompat.sortedToc())
    {
        const auto& iter = *argList::ignoreOptionsCompat.cfind(k);

        const bool hasArg = iter.first;
        const int until = abs(iter.second);

        Info<< ".TP\n\\fB\\-" << k << "\\fR";

        if (hasArg)
        {
            Info<< " \\fIarg\\fR";
        }

        Info<< nl << "This option is ignored";

        if (until)
        {
            Info<< " after " << until << ".";
        }
        Info<< nl;
    }
}


void Foam::argList::printCompat() const
{
    const label nopt
    (
        argList::validOptionsCompat.size()
      + argList::ignoreOptionsCompat.size()
    );

    Info<< nopt << " compatibility options for " << executable_ << nl;

    if (!nopt)
    {
        return;
    }

    const int col1(32), col2(32);

    Info<< nl
        << "|" << setf(ios_base::left) << setw(col1) << " Old option"
        << "|" << setf(ios_base::left) << setw(col2) << " New option"
        << "| Comment" << nl;

    Info<< setfill('-');
    Info<< "|" << setf(ios_base::left) << setw(col1) << ""
        << "|" << setf(ios_base::left) << setw(col2) << ""
        << "|------------" << nl;

    Info<< setfill(token::SPACE);

    for (const word& k : argList::validOptionsCompat.sortedToc())
    {
        const auto& iter = *argList::validOptionsCompat.cfind(k);

        const word& optName = iter.first;
        const int until = abs(iter.second);

        Info<< "| -" << setf(ios_base::left) << setw(col1-2) << k
            << "| -" << setf(ios_base::left) << setw(col2-2) << optName
            << "|";

        if (until)
        {
            Info<< " until " << until;
        }
        Info<< nl;
    }

    for (const word& k : argList::ignoreOptionsCompat.sortedToc())
    {
        const auto& iter = *argList::ignoreOptionsCompat.cfind(k);

        const bool hasArg = iter.first;
        const int until = abs(iter.second);

        Info<< "| -" << setf(ios_base::left) << setw(col1-2);

        if (hasArg)
        {
            Info<< (k + " <arg>").c_str();
        }
        else
        {
            Info<< k;
        }

        Info<< "| ";
        Info<< setf(ios_base::left) << setw(col2-1) << "ignored" << "|";
        if (until)
        {
            Info<< " after " << until;
        }
        Info<< nl;
    }

    Info<< setfill('-');
    Info<< "|" << setf(ios_base::left) << setw(col1) << ""
        << "|" << setf(ios_base::left) << setw(col2) << ""
        << "|------------" << nl;

    Info<< setfill(token::SPACE);
}


// ************************************************************************* //

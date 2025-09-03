/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "entry.H"
#include "dictionary.H"
#include "SpanStream.H"
#include "JobInfo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::entry::disableFunctionEntries
(
    Foam::debug::infoSwitch("disableFunctionEntries", 0)
);


Foam::entry::inputMode Foam::entry::globalInputMode = inputMode::MERGE;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::entry::reportReadWarning
(
    const IOstream& is,
    const std::string& msg
)
{
    std::cerr
        << "--> FOAM Warning :\n"
        << "    Reading \"" << is.relativeName()
        << "\" at line " << is.lineNumber() << '\n'
        << "    " << msg << std::endl;
}


void Foam::entry::resetInputMode()
{
    globalInputMode = inputMode::MERGE;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entry::entry(const keyType& keyword)
:
    IDLList<entry>::link(),
    keyword_(keyword)
{}


Foam::entry::entry(const entry& e)
:
    IDLList<entry>::link(),
    keyword_(e.keyword_)
{}


Foam::autoPtr<Foam::entry> Foam::entry::clone() const
{
    return clone(dictionary::null);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::entry::raiseBadInput(const ITstream& is) const
{
    const word& keyword = keyword_;

    // Can use FatalIOError instead of SafeFatalIOError
    // since predicate checks are not used at the earliest stages
    FatalIOError
    (
        "",                 // functionName
        "",                 // sourceFileName
        0,                  // sourceFileLineNumber
        this->relativeName(), // ioFileName
        is.lineNumber()     // ioStartLineNumber
    )
        << "Entry '" << keyword << "' with invalid input" << nl << nl
        << exit(FatalIOError);
}


void Foam::entry::checkITstream(const ITstream& is) const
{
    const label remaining = (is.size() ? is.nRemainingTokens() : -100);

    if (!remaining)
    {
        return;
    }

    const word& keyword = keyword_;

    // Similar to SafeFatalIOError
    if (JobInfo::constructed)
    {
        OSstream& err =
            FatalIOError
            (
                "",                 // functionName
                "",                 // sourceFileName
                0,                  // sourceFileLineNumber
                this->relativeName(), // ioFileName
                is.lineNumber()     // ioStartLineNumber
            );

        if (remaining > 0)
        {
            err
                << "Entry '" << keyword
                << "' has " << remaining << " excess tokens in stream"
                << nl << nl
                << "    ";
            is.writeList(err, 0);  // <- flatOutput
        }
        else
        {
            err
                << "Entry '" << keyword
                << "' had no tokens in stream"
                << nl << nl;
        }

        err << exit(FatalIOError);
    }
    else
    {
        // Not yet constructed

        std::cerr
            << nl
            << "--> FOAM FATAL IO ERROR:" << nl;

        if (remaining > 0)
        {
            std::cerr
                << "Entry '" << keyword << "' has "
                << remaining << " excess tokens in stream" << nl << nl;
        }
        else
        {
            std::cerr
                << "Entry '" << keyword
                << "' had no tokens in stream" << nl << nl;
        }

        std::cerr
            << "file: " << this->relativeName()
            << " at line " << is.lineNumber() << '.' << nl
            << std::endl;

        std::exit(1);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::entry::operator=(const entry& e)
{
    if (this == &e)
    {
        return;  // Self-assignment is a no-op
    }

    keyword_ = e.keyword_;
}


bool Foam::entry::operator==(const entry& e) const
{
    if (this == &e)
    {
        return true;
    }
    if (keyword_ != e.keyword_)
    {
        return false;
    }

    // Compare contents (as strings)
    OCharStream content1;
    content1 << *this;

    OCharStream content2;
    content2 << e;

    return
    (
        content1.view().size() == content2.view().size()
     &&
        std::equal
        (
            content1.view().begin(),
            content1.view().end(),
            content2.view().begin()
        )
    );
}


bool Foam::entry::operator!=(const entry& e) const
{
    return !operator==(e);
}


// ************************************************************************* //

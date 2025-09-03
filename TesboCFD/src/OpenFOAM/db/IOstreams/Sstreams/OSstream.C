/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "token.H"
#include "OSstream.H"
#include <algorithm>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::OSstream::write(const token& tok)
{
    // Direct token handling only for some types

    switch (tok.type())
    {
        case token::tokenType::FLAG :
        {
            // silently consume the flag
            return true;
        }

        case token::tokenType::DIRECTIVE :
        {
            // Token stored with leading '#' sigil - output directly
            write(tok.wordToken());
            return true;
        }

        case token::tokenType::EXPRESSION :
        {
            // Token stored with surrounding '${{ .. }}' - output directly
            writeQuoted(tok.stringToken(), false);
            return true;
        }

        case token::tokenType::VARIABLE :
        {
            // Token stored with leading '$' sigil - output directly
            writeQuoted(tok.stringToken(), false);
            return true;
        }

        case token::tokenType::VERBATIM :
        {
            // Token stored without surrounding '#{ .. #}'. Add on output
            write(char(token::HASH));
            write(char(token::BEGIN_BLOCK));
            writeQuoted(tok.stringToken(), false);
            write(char(token::HASH));
            write(char(token::END_BLOCK));

            return true;
        }

        case token::tokenType::CHAR_DATA :
        {
            // Character content - written without quotes
            writeQuoted(tok.stringToken(), false);
            return true;
        }

        default:
            break;
    }

    return false;
}


Foam::Ostream& Foam::OSstream::write(const char c)
{
    os_ << c;
    syncState();

    // Advance line number on newline
    if (c == token::NL) ++lineNumber_;
    return *this;
}


Foam::Ostream& Foam::OSstream::writeQuoted
(
    const char* str,
    std::streamsize len,
    const bool quoted
)
{
    if (!str || len <= 0) return *this;

    const char* last = (str + len);

    if (!quoted)
    {
        #if __cplusplus >= 201703L
        os_ << std::string_view(str, len);
        #else
        for (const char* iter = str; iter != last; ++iter)
        {
            os_ << *iter;
        }
        #endif
        syncState();

        // Unquoted, only advance line number on newline
        lineNumber_ += std::count(str, last, '\n');
        return *this;
    }


    // Output with surrounding quotes and backslash escaping
    // - functionality like std::quoted (from <iomanip>), while also
    //   counting the newlines.

    os_ << token::DQUOTE;

    unsigned backslash = 0;
    for (auto iter = str; iter != last; ++iter)
    {
        const char c = *iter;

        if (c == '\\')
        {
            ++backslash;
            continue;       // Delay output until escaped character is known
        }
        else if (c == token::NL)
        {
            ++backslash;    // Add backslash escape for newline
            ++lineNumber_;  // Advance line number on newline
        }
        else if (c == token::DQUOTE)
        {
            ++backslash;    // Add backslash escape for quote
        }

        // output all pending backslashes
        while (backslash)
        {
            os_ << '\\';
            --backslash;
        }

        os_ << c;
    }

    // silently drop any trailing backslashes
    // they would otherwise appear like an escaped end-quote
    os_ << token::DQUOTE;

    syncState();
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const char* str)
{
    if (!str) return *this;

    const char* last = (str + strlen(str));

    os_ << str;
    syncState();

    // Advance line number on newline
    lineNumber_ += std::count(str, last, '\n');
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const word& str)
{
    // Unquoted, and no newlines expected.
    os_ << str;
    syncState();

    return *this;
}


Foam::Ostream& Foam::OSstream::write(const std::string& str)
{
    return writeQuoted(str.data(), str.size(), true);
}


Foam::Ostream& Foam::OSstream::write(const int32_t val)
{
    os_ << val;
    syncState();
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const int64_t val)
{
    os_ << val;
    syncState();
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const float val)
{
    os_ << val;
    syncState();
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const double val)
{
    os_ << val;
    syncState();
    return *this;
}


Foam::Ostream& Foam::OSstream::write(const char* data, std::streamsize count)
{
    beginRawWrite(count);
    writeRaw(data, count);
    endRawWrite();

    return *this;
}


bool Foam::OSstream::beginRawWrite(std::streamsize count)
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalIOErrorInFunction(*this)
            << "stream format not binary"
            << abort(FatalIOError);
    }

    os_ << token::BEGIN_LIST;
    syncState();
    return os_.good();
}


bool Foam::OSstream::endRawWrite()
{
    os_ << token::END_LIST;
    syncState();
    return os_.good();
}


Foam::Ostream& Foam::OSstream::writeRaw
(
    const char* data,
    std::streamsize count
)
{
    // No check for IOstreamOption::BINARY since this is either done in the
    // beginRawWrite() method, or the caller knows what they are doing.

    os_.write(data, count);
    syncState();
    return *this;
}


void Foam::OSstream::indent()
{
    for (unsigned short i = 0; i < indentLevel_*indentSize_; ++i)
    {
        os_ << ' ';
    }
    syncState();
}


void Foam::OSstream::flush()
{
    os_.flush();
}


void Foam::OSstream::endl()
{
    write('\n');
    os_.flush();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char Foam::OSstream::fill() const
{
    return os_.fill();
}


char Foam::OSstream::fill(const char fillch)
{
    return os_.fill(fillch);
}


int Foam::OSstream::width() const
{
    return os_.width();
}


int Foam::OSstream::width(const int w)
{
    return os_.width(w);
}


int Foam::OSstream::precision() const
{
    return os_.precision();
}


int Foam::OSstream::precision(const int p)
{
    return os_.precision(p);
}


// ************************************************************************* //

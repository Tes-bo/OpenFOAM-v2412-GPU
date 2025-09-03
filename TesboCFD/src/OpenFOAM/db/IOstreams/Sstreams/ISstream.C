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

#include "ISstream.H"
#include "int.H"
#include "token.H"
#include <cctype>
#include <cstring>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Truncate error message for readability
static constexpr const unsigned errLen = 80;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Convert a single character to a word with length 1
inline Foam::word charToWord(char c)
{
    return Foam::word(std::string(1, c), false);
}


// Permit slash-scoping of entries
inline bool validVariableChar(char c)
{
    return (Foam::word::valid(c) || c == '/');
}


inline void inplaceTrimRight(std::string& s)
{
    auto end = s.length();
    if (end)
    {
        while (end && Foam::isspace(s[end-1]))
        {
            --end;
        }

        s.erase(end);
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::ISstream::readCompoundToken
(
    token& tok,
    const word& compoundType
)
{
    return tok.readCompoundToken(compoundType, *this);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ISstream::seekCommentEnd_Cstyle()
{
    // Search for end of C-style comment - "*/"

    // Can use getLine(nullptr, '*') in the logic,
    // but written out looks less obscure

    char c = 0;
    bool star = false;

    while (get(c))
    {
        if (c == '*')
        {
            star = true;
        }
        else if (star)
        {
            star = false;
            if (c == '/')
            {
                // Matched "*/"
                return true;
            }
        }
    }

    // Exhausted stream without finding "*/" sequence
    return false;
}


char Foam::ISstream::nextValid()
{
    char c = 0;

    // Get next non-whitespace character
    while (get(c))
    {
        if (isspace(c))
        {
            continue;
        }

        // Check if this starts a C/C++ comment
        if (c == '/')
        {
            if (!get(c))
            {
                // Cannot get another character - return this one
                return '/';
            }

            if (c == '/')
            {
                // C++ comment: discard through newline
                (void) getLine(nullptr, '\n');
            }
            else if (c == '*')
            {
                // C-style comment: discard through to "*/" ending
                if (!seekCommentEnd_Cstyle())
                {
                    return 0;  // Premature end of stream
                }
            }
            else
            {
                // The '/' did not start a C/C++ comment - return it
                putback(c);
                return '/';
            }
        }
        else
        {
            // A valid character - return it
            return c;
        }
    }

    return 0;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Read a verbatim string (excluding block delimiters),
// continuing until a closing "#}" has been found.
//
// The leading "#{" removed from stream prior to calling.
static ISstream& readVerbatim
(
    ISstream& is,
    std::string& str
)
{
    constexpr const unsigned bufLen = 8000;
    static char buf[bufLen];

    unsigned nChar = 0;
    char c;

    str.clear();
    while (is.get(c))
    {
        if (c == token::HASH)
        {
            char nextC;
            is.get(nextC);
            if (nextC == token::END_BLOCK)
            {
                // Found closing "#}" sequence
                str.append(buf, nChar);
                return is;
            }
            else
            {
                // Re-analyze the character
                is.putback(nextC);
            }
        }

        buf[nChar++] = c;
        if (nChar == bufLen)  // Flush full buffer
        {
            str.append(buf, nChar);
            nChar = 0;
        }
    }


    // Abnormal exit of the loop
    str.append(buf, nChar);  // Finalize pending content
    strncpy(buf, str.c_str(), errLen);
    buf[errLen] = '\0';

    FatalIOErrorInFunction(is)
        << "Problem while reading verbatim \"" << buf
        << "...\" [after " << str.length() << " chars]\n"
        << exit(FatalIOError);

    return is;
}


// Read a variable or expression.
// Handles "$var" and "${var}" forms, permits '/' scoping character.
// Also handles "${{expr}}".
//
// Return the token type or ERROR
//
// The leading "${" or "$c" removed from stream prior to calling.
static token::tokenType readVariable
(
    ISstream& is,
    std::string& str,
    char c  // Next character after '$'
)
{
    constexpr const unsigned bufLen = 1024;
    static char buf[bufLen];

    token::tokenType tokType(token::tokenType::VARIABLE);

    // The first two characters are known:
    buf[0] = token::DOLLAR;
    buf[1] = c;

    unsigned nChar = 2; // Starts with two characters
    unsigned depth = 0; // Depth of {..} nesting

    str.clear();
    if (c == token::BEGIN_BLOCK)
    {
        // Processing '${variable}' or '${{expr}}'
        ++depth;

        int lookahead = is.peek();
        if (lookahead == token::BEGIN_BLOCK)
        {
            // Looks like '${{expr...'
            tokType = token::tokenType::EXPRESSION;
        }
        else if (lookahead == token::END_BLOCK)
        {
            // Looks like '${}'
            IOWarningInFunction(is)
                << "Ignoring empty ${}" << endl;
            return token::tokenType::ERROR;
        }

        while (is.get(c))
        {
            buf[nChar++] = c;

            if (c == token::BEGIN_BLOCK)
            {
                ++depth;
            }
            else if (c == token::END_BLOCK)
            {
                --depth;
                if (!depth)
                {
                    // Found closing '}' character
                    str.append(buf, nChar);
                    return tokType;
                }
            }
            else if (c == '/' && tokType == token::tokenType::EXPRESSION)
            {
                // Strip C/C++ comments from expressions
                // Note: could also peek instead of get/putback

                if (!is.get(c))
                {
                    break;  // Premature end of stream
                }
                else if (c == '/')
                {
                    --nChar;  // Remove initial '/' from buffer

                    // C++ comment: discard through newline
                    (void) is.getLine(nullptr, '\n');
                }
                else if (c == '*')
                {
                    --nChar;  // Remove initial '/' from buffer

                    // C-style comment: seek "*/" ending
                    if (!is.seekCommentEnd_Cstyle())
                    {
                        break;  // Premature end of stream
                    }
                }
                else
                {
                    // Re-analyze the character
                    is.putback(c);
                }
            }

            if (nChar == bufLen)  // Flush full buffer
            {
                str.append(buf, nChar);
                nChar = 0;
            }
        }


        // Abnormal exit of the loop

        str.append(buf, nChar);  // Finalize pending content
        strncpy(buf, str.c_str(), errLen);
        buf[errLen] = '\0';

        FatalIOErrorInFunction(is)
            << "stream terminated while reading variable '" << buf
            << "...' [after " << str.length() << " chars]\n"
            << exit(FatalIOError);

        return token::tokenType::ERROR;
    }
    else if (validVariableChar(c))
    {
        // Processing '$variable'

        while (is.get(c))
        {
            if (!validVariableChar(c))
            {
                is.putback(c);
                break;
            }

            if (c == token::BEGIN_LIST)
            {
                ++depth;
            }
            else if (c == token::END_LIST)
            {
                if (!depth)
                {
                    // Closed ')' without opening '(':
                    // - don't consider it part of our input
                    is.putback(c);
                    break;
                }
                --depth;
            }

            buf[nChar++] = c;
            if (nChar == bufLen)  // Flush full buffer
            {
                str.append(buf, nChar);
                nChar = 0;
            }
        }

        str.append(buf, nChar);  // Finalize pending content

        if (depth)
        {
            strncpy(buf, str.c_str(), errLen);
            buf[errLen] = '\0';

            IOWarningInFunction(is)
                << "Missing " << depth
                << " closing ')' while parsing" << nl << nl
                << buf << endl;
        }

        return tokType;
    }
    else
    {
        // Invalid character. Terminate string (for message)

        buf[nChar--] = '\0';

        IOWarningInFunction(is)
            << "Ignoring bad variable name: " << buf << nl << endl;
    }

    return token::tokenType::ERROR;
}


// Raw, low-level get into a string.
// Continues reading after an initial opening delimiter (eg, '{')
// until it finds the matching closing delimiter (eg, '}')
static bool readUntilBalancedDelimiter
(
    ISstream& is,
    std::string& str,
    const bool stripComments,
    const char delimOpen,
    const char delimClose
)
{
    constexpr const unsigned bufLen = 1024;
    static char buf[bufLen];

    unsigned nChar = 0;
    unsigned depth = 1;  // Initial '{' already seen by caller
    char c = 0;

    str.clear();
    while (is.get(c))
    {
        if ((str.empty() && !nChar) && isspace(c))
        {
            continue;  // Ignore leading whitespace
        }

        buf[nChar++] = c;

        // Note: no '\' escape handling needed at the moment

        if (c == delimOpen)
        {
            ++depth;
        }
        else if (c == delimClose)
        {
            --depth;
            if (!depth)
            {
                // Closing character - do not include in output
                --nChar;
                str.append(buf, nChar);
                inplaceTrimRight(str);  // Remove trailing whitespace
                return true;
            }
        }
        else if (stripComments && c == '/')
        {
            // Strip C/C++ comments from expressions
            // Note: could also peek instead of get/putback

            if (!is.get(c))
            {
                break;  // Premature end of stream
            }
            else if (c == '/')
            {
                --nChar;  // Remove initial '/' from buffer

                // C++ comment: discard through newline
                (void) is.getLine(nullptr, '\n');
            }
            else if (c == '*')
            {
                --nChar;  // Remove initial '/' from buffer

                // C-style comment: discard through to "*/" ending
                if (!is.seekCommentEnd_Cstyle())
                {
                    break;  // Premature end of stream
                }
            }
            else
            {
                // Reanalyze the char
                is.putback(c);
            }
        }

        if (nChar == bufLen)
        {
            str.append(buf, nChar);  // Flush full buffer
            nChar = 0;
        }
    }


    // Abnormal exit of the loop

    str.append(buf, nChar);  // Finalize pending content
    inplaceTrimRight(str);   // Remove trailing whitespace

    // Exhausted stream without finding closing sequence
    return false;
}

} // End namespace Foam


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::ISstream::continueReadUntilRightBrace
(
    std::string& str,
    const bool stripComments
)
{
    return
        readUntilBalancedDelimiter
        (
            *this,
            str,
            stripComments,
            token::BEGIN_BLOCK,
            token::END_BLOCK
        );
}


Foam::Istream& Foam::ISstream::read(token& t)
{
    constexpr const unsigned bufLen = 128; // Max length for labels/scalars
    static char buf[bufLen];

    // Return the putback token if it exists
    if (Istream::getBack(t))
    {
        return *this;
    }

    // Reset token, adjust its line number according to the stream
    t.reset();
    t.lineNumber(this->lineNumber());

    // Assume that the streams supplied are in working order.
    // Lines are counted by '\n'

    // Get next 'valid character': i.e. proceed through any whitespace
    // and/or comments until a semantically valid character is found

    char c = nextValid();

    // Return on error
    if (!c)
    {
        t.setBad();
        return *this;
    }

    // Analyse input starting with this character.
    switch (c)
    {
        // Check for punctuation first - same as token::isseparator()

        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::PLUS :
        // NB: token::MINUS handled later as the possible start of a Number
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // String: enclosed by double quotes.
        case token::DQUOTE :
        {
            putback(c);

            string val;
            if (read(val).bad())
            {
                t.setBad();
            }
            else
            {
                t = std::move(val); // Move contents to token
            }

            return *this;
        }

        // Verbatim string '#{ .. #}' or dictionary '#directive'
        case token::HASH :
        {
            char nextC;
            int lookahead = peek();

            if (lookahead == token::BEGIN_BLOCK)
            {
                // Verbatim string: #{ ... #}
                // Token stored without the surrounding delimiters

                (void) get(nextC);  // Discard '{' lookahead

                string val;
                if (readVerbatim(*this, val).bad())
                {
                    t.setBad();
                }
                else
                {
                    t = std::move(val); // Move contents to token
                    t.setType(token::tokenType::VERBATIM);
                }
            }
            else if (read(nextC).bad())
            {
                // Return lone '#' as word
                t = charToWord(c);
            }
            else if (word::valid(nextC))
            {
                // Directive (wordToken) beginning with '#'. Eg, "#include"
                // Put back both so that '#...' is included in the directive

                putback(nextC);
                putback(c);

                word val;
                if (read(val).bad())
                {
                    t.setBad();
                }
                else
                {
                    t = std::move(val); // Move contents to token
                    t.setType(token::tokenType::DIRECTIVE);
                }
            }
            else
            {
                // '#' followed by non-word. Just ignore leading '#'?
                putback(nextC);

                IOWarningInFunction(*this)
                    << "Invalid sequence #" << char(nextC)
                    << " ... ignoring the leading '#'" << nl << endl;
            }

            return *this;
        }

        // Dictionary variable or ${{ expression }}
        case token::DOLLAR :
        {
            char nextC;
            if (read(nextC).bad())
            {
                // Return lone '$' as word. Could also ignore
                t = charToWord(c);
            }
            else
            {
                // NB: the parser is slightly generous here.
                // It will also accept '$  {' as input.
                // - to be revisited (2021-05-17)

                string val;
                token::tokenType tokType = readVariable(*this, val, nextC);
                if (tokType == token::tokenType::ERROR)
                {
                    t.setBad();
                }
                else
                {
                    t = std::move(val); // Move contents to token
                    t.setType(tokType);
                }
            }

            return *this;
        }

        // Number: integer or floating point
        //
        // ideally match the equivalent of this regular expression
        //
        //    /[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([Ee][-+]?[0-9]+)?/
        //
        case '-' :
        case '.' :
        case '0' : case '1' : case '2' : case '3' : case '4' :
        case '5' : case '6' : case '7' : case '8' : case '9' :
        {
            label labelVal = (c != '.'); // used as bool here

            unsigned nChar = 0;
            buf[nChar++] = c;

            // get everything that could resemble a number and let
            // readScalar determine the validity
            while
            (
                is_.get(c)
             && (
                    isdigit(c)
                 || c == '+'
                 || c == '-'
                 || c == '.'
                 || c == 'E'
                 || c == 'e'
                )
            )
            {
                if (labelVal)
                {
                    labelVal = isdigit(c);
                }

                buf[nChar++] = c;
                if (nChar == bufLen)
                {
                    // Runaway argument - avoid buffer overflow
                    buf[bufLen-1] = '\0';

                    FatalIOErrorInFunction(*this)
                        << "Number '" << buf << "...'\n"
                        << "    is too long (max. " << bufLen << " characters)"
                        << exit(FatalIOError);

                    t.setBad();
                    return *this;
                }
            }
            buf[nChar] = '\0';  // Terminate string

            syncState();

            if (is_.bad())
            {
                t.setBad();
            }
            else
            {
                is_.putback(c);

                if (nChar == 1 && buf[0] == '-')
                {
                    // A single '-' is punctuation
                    t = token::punctuationToken(token::MINUS);
                }
                else if (labelVal && Foam::read(buf, labelVal))
                {
                    t = labelVal;
                }
                else
                {
                    scalar scalarVal;

                    if (readScalar(buf, scalarVal))
                    {
                        // A scalar or too big to fit as a label
                        t = scalarVal;
                    }
                    else
                    {
                        t.setBad();
                    }
                }
            }

            return *this;
        }

        // Should be a word (which can also be a single character)
        default:
        {
            putback(c);

            word val;
            if (read(val).bad())
            {
                t.setBad();
            }
            else if
            (
                !token::compound::isCompound(val)
             || !readCompoundToken(t, val)
            )
            {
                t = std::move(val);  // Move contents to token
            }

            return *this;
        }
    }
}


Foam::Istream& Foam::ISstream::read(char& c)
{
    c = nextValid();
    return *this;
}


Foam::Istream& Foam::ISstream::read(word& str)
{
    constexpr const unsigned bufLen = 1024;
    static char buf[bufLen];

    unsigned nChar = 0;
    unsigned depth = 0;  // Depth of (..) nesting
    char c;

    str.clear();
    while (get(c))
    {
        if (!word::valid(c))
        {
            putback(c);
            break;
        }

        if (c == token::BEGIN_LIST)
        {
            ++depth;
        }
        else if (c == token::END_LIST)
        {
            if (!depth)
            {
                // Closed ')' without opening '(':
                // - don't consider it part of our input
                putback(c);
                break;
            }
            --depth;
        }

        buf[nChar++] = c;
        if (nChar == bufLen)  // Flush full buffer
        {
            str.append(buf, nChar);
            nChar = 0;
        }
    }

    str.append(buf, nChar);  // Finalize pending content

    if (bad())
    {
        // Could probably skip this check

        strncpy(buf, str.c_str(), errLen);
        buf[errLen] = '\0';

        FatalIOErrorInFunction(*this)
            << "Problem while reading word '" << buf
            << "...' [after " << str.length() << " chars]\n"
            << exit(FatalIOError);

        return *this;
    }

    if (str.empty())
    {
        FatalIOErrorInFunction(*this)
            << "Invalid first character found : " << c
            << exit(FatalIOError);
    }
    else if (depth)
    {
        strncpy(buf, str.c_str(), errLen);
        buf[errLen] = '\0';

        IOWarningInFunction(*this)
            << "Missing " << depth
            << " closing ')' while parsing" << nl << nl
            << buf << nl << endl;
    }

    return *this;
}


Foam::Istream& Foam::ISstream::read(string& str)
{
    constexpr const unsigned bufLen = 1024;
    static char buf[bufLen];

    unsigned nChar = 0;
    char c;

    if (!get(c))
    {
        FatalIOErrorInFunction(*this)
            << "cannot read start of string"
            << exit(FatalIOError);

        return *this;
    }

    // Note, we could also handle single-quoted strings here (if desired)
    if (c != token::DQUOTE)
    {
        FatalIOErrorInFunction(*this)
            << "Incorrect start of string character found : " << c
            << exit(FatalIOError);

        return *this;
    }

    str.clear();
    bool escaped = false;
    while (get(c))
    {
        if (c == '\\')
        {
            escaped = !escaped;  // Toggle state (retains backslashes)
        }
        else if (c == token::DQUOTE)
        {
            if (escaped)
            {
                escaped = false;
                --nChar;  // Overwrite backslash
            }
            else
            {
                // Done reading
                str.append(buf, nChar);
                return *this;
            }
        }
        else if (c == token::NL)
        {
            if (escaped)
            {
                escaped = false;
                --nChar;  // Overwrite backslash
            }
            else
            {
                str.append(buf, nChar);  // Finalize pending content
                strncpy(buf, str.c_str(), errLen);
                buf[errLen] = '\0';

                FatalIOErrorInFunction(*this)
                    << "Unescaped '\\n' while reading string \"" << buf
                    << "...\" [after " << str.length() << " chars]\n"
                    << exit(FatalIOError);

                return *this;
            }
        }
        else
        {
            escaped = false;
        }

        buf[nChar++] = c;
        if (nChar == bufLen)  // Flush full buffer
        {
            // Keep lookback character (eg, for backslash escaping)
            str.append(buf, nChar-1);
            nChar = 1;
            buf[0] = c;
        }
    }


    // Abnormal exit of the loop
    // Don't worry about a dangling backslash if string terminated prematurely

    str.append(buf, nChar);  // Finalize pending content
    strncpy(buf, str.c_str(), errLen);
    buf[errLen] = '\0';

    FatalIOErrorInFunction(*this)
        << "Problem while reading string \"" << buf << "...\""
        << exit(FatalIOError);

    return *this;
}


Foam::Istream& Foam::ISstream::read(label& val)
{
    is_ >> val;
    syncState();
    return *this;
}


Foam::Istream& Foam::ISstream::read(float& val)
{
    is_ >> val;
    syncState();
    return *this;
}


Foam::Istream& Foam::ISstream::read(double& val)
{
    is_ >> val;
    syncState();
    return *this;
}


Foam::Istream& Foam::ISstream::read(char* data, std::streamsize count)
{
    beginRawRead();
    readRaw(data, count);
    endRawRead();

    return *this;
}


Foam::Istream& Foam::ISstream::readRaw(char* data, std::streamsize count)
{
    if (count)
    {
        if (data)
        {
            is_.read(data, count);
        }
        else
        {
            is_.ignore(count);
        }
    }
    syncState();
    return *this;
}


bool Foam::ISstream::beginRawRead()
{
    if (format() != IOstreamOption::BINARY)
    {
        FatalIOErrorInFunction(*this)
            << "stream format not binary"
            << exit(FatalIOError);
    }

    readBegin("binaryBlock");
    syncState();
    return is_.good();
}


bool Foam::ISstream::endRawRead()
{
    readEnd("binaryBlock");
    syncState();
    return is_.good();
}


void Foam::ISstream::rewind()
{
    lineNumber_ = 1;      // Reset line number

    stdStream().clear();  // Clear the iostate error state flags
    setGood();            // Sync local copy of iostate

    stdStream().rdbuf()->pubseekpos(0, std::ios_base::in);

    // NOTE: this form of rewind does not work with igzstream.
    // However, igzstream is usually wrapped as IFstream which has its
    // own dedicated rewind treatment for igzstream.
}


// ************************************************************************* //

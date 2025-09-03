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

#include "stringOps.H"
#include "typeInfo.H"
#include "etcFiles.H"
#include "UPstream.H"
#include "StringStream.H"
#include "OSstream.H"
#include "OSspecific.H"
#include <cctype>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Return the file location mode (string) as a numerical value.
//
// - u : location mask 0700
// - g : location mask 0070
// - o : location mask 0007
// - a : location mask 0777
//
static inline unsigned short modeToLocation
(
    const std::string& mode,
    std::size_t pos = 0
)
{
    unsigned short where(0);

    if (std::string::npos != mode.find('u', pos)) { where |= 0700; } // User
    if (std::string::npos != mode.find('g', pos)) { where |= 0070; } // Group
    if (std::string::npos != mode.find('o', pos)) { where |= 0007; } // Other
    if (std::string::npos != mode.find('a', pos)) { where |= 0777; } // All

    return where;
}


// Expand a leading <tag>/
// Convenient for frequently used directories
//
//   <etc>/        => user/group/other etc - findEtcEntry()
//   <etc(:[ugoa]+)?>/ => user/group/other etc - findEtcEntry()
//   <case>/       => FOAM_CASE directory
//   <constant>/   => FOAM_CASE/constant directory
//   <system>/     => FOAM_CASE/system directory
static void expandLeadingTag(std::string& s, const char b, const char e)
{
    if (s[0] != b)
    {
        return;
    }

    auto delim = s.find(e);
    if (std::string::npos == delim)
    {
        return;  // Error: no closing delim - ignore expansion
    }

    fileName file;

    const char nextC = s[++delim];

    // Require the following character to be '/' or the end of string.
    if (nextC)
    {
        if (nextC != '/')
        {
            return;
        }

        file.assign(s.substr(delim + 1));
    }

    const std::string tag(s, 1, delim-2);
    const auto tagLen = tag.length();

    // Note that file is also allowed to be an empty string.

    if (tag == "etc")
    {
        s = findEtcEntry(file);
    }
    else if (tag == "case")
    {
        s = fileName(Foam::getEnv("FOAM_CASE"))/file;
    }
    else if (tag == "constant" || tag == "system")
    {
        s = fileName(Foam::getEnv("FOAM_CASE"))/tag/file;
    }
    else if (tagLen >= 4 && tag.compare(0, 4, "etc:") == 0)
    {
        // <etc:[ugoa]+> type of tag - convert "ugo" to numeric

        s = findEtcEntry(file, modeToLocation(tag,4));
    }
}


// Expand a leading tilde
//   ~/        => home directory
//   ~user     => home directory for specified user
//   Deprecated ~OpenFOAM => <etc> instead
static void expandLeadingTilde(std::string& s)
{
    if (s[0] != '~')
    {
        return;
    }

    std::string user;
    fileName file;

    const auto slash = s.find('/');
    if (slash == std::string::npos)
    {
        user = s.substr(1);
    }
    else
    {
        user = s.substr(1, slash - 1);
        file = s.substr(slash + 1);
    }

    // NB: be a bit lazy and expand ~unknownUser as an
    // empty string rather than leaving it untouched.
    // otherwise add extra test

    if (user == "OpenFOAM")
    {
        // Compat Warning
        const int version(1806);

        if (error::master())
        {
            std::cerr
                << nl
                << "--> FOAM Warning :" << nl
                << "    Found [v" << version << "] '"
                << "~OpenFOAM" << "' string expansion instead of '"
                << "<etc>" << "' in string\n\"" << s << "\"\n" << nl
                << std::endl;

            error::warnAboutAge("expansion", version);
        }

        s = findEtcFile(file);
    }
    else
    {
        s = home(user)/file;
    }
}


// Expand leading contents:  "./", "~..", "<tag>/"
static void expandLeading(std::string& s)
{
    if (s.empty())
    {
        return;
    }

    switch (s[0])
    {
        case '.':
        {
            // Expand a lone '.' and an initial './' into cwd
            if (s.size() == 1)
            {
                s = cwd();
            }
            else if (s[1] == '/')
            {
                s.replace(0, 1, cwd());
            }
            break;
        }
        case '<':
        {
            expandLeadingTag(s, '<', '>');
            break;
        }
        case '~':
        {
            expandLeadingTilde(s);
            break;
        }
    }
}


// Serialize an entry (primitive or dictionary) with special treatment
// for primitive entries that are already a string-type.
static inline std::string entryToString
(
    const entry* eptr,
    const bool allowSubDict
)
{
    std::string str;

    if (eptr)
    {
        // Note, for OpenFOAM-v2212 and earlier: used 'fixed' notation
        // to force floating point numbers to be printed with at least
        // some decimal digits. However, in the meantime we are more
        // flexible with handling float/int input so remove this constraint.

        if (eptr->isDict())
        {
            if (allowSubDict)
            {
                OStringStream buf;
                buf.precision(16);  // Some reasonably high precision

                eptr->dict().write(buf, false);
                str = buf.str();
            }
            else
            {
                // Ignore silently...
            }
        }
        else
        {
            // Serialized with spaces (primitiveEntry)
            ITstream& its = eptr->stream();
            return its.toString();
        }
    }

    return str;
}

} // End namespace Foam


// Details for handling dictionary expansion

namespace
{

// Acceptable values for $variable names.
//
// Similar to word::valid(), except we don't have the benefit of a parser
// to filter out other unacceptable entries for us.
//
// Does not currently accept '/' in a variable name.
// We would like "$file/$name" to expand as two variables.
static inline bool validVariableChar(char c)
{
    return
    (
        std::isalnum(c)
     || c == '.'
     || c == ':'
     || c == '_'
    );
}


//  Find the type/position of the ":-" or ":+" alternative values
//  Returns 0, '-', '+' corresponding to not-found or ':-' or ':+'
static inline int findParameterAlternative
(
    const std::string& s,
    std::string::size_type& pos,
    std::string::size_type endPos
)
{
    while (pos != std::string::npos)
    {
        pos = s.find(':', pos);
        if (pos != std::string::npos)
        {
            if (pos < endPos)
            {
                // in-range: check for '+' or '-' following the ':'
                const int altType = s[pos+1];
                if (altType == '+' || altType == '-')
                {
                    return altType;
                }

                ++pos;    // unknown/unsupported - continue at next position
            }
            else
            {
                // out-of-range: abort
                pos = std::string::npos;
            }
        }
    }

    return 0;
}


// For input string of "$variable with other" return the length of
// the variable.
//
// Intentionally will not capture ':+', ':-' alterations. Use ${ .. } for that
static inline std::string::size_type findVariableLen
(
    const std::string& s,
    std::string::size_type pos,
    const char sigil = '$'
)
{
    std::string::size_type len = 0;

    if (pos < s.length())
    {
        if (s[pos] == sigil)
        {
            // Skip leading '$' in the count!
            ++pos;
        }

        for
        (
            auto iter = s.cbegin() + pos;
            iter != s.cend() && validVariableChar(*iter);
            ++iter
        )
        {
            ++len;
        }
    }

    return len;
}

} // End namespace anonymous


namespace Foam
{

// Get dictionary or (optionally) environment variable
//
// Handles default and alternative values as per the POSIX shell.
//  \code
//      ${parameter:-defValue}
//      ${parameter:+altValue}
//  \endcode
static Foam::string getVariable
(
    const word& name,
    const dictionary* dictptr,
    const bool allowEnv,
    const bool allowEmpty,
    const bool allowSubDict
)
{
    // The type/position of the ":-" or ":+" alternative values
    std::string::size_type altPos = 0;

    // Check for parameter:-word or parameter:+word
    const int altType =
        findParameterAlternative(name, altPos, name.size()-1);

    const word lookupName =
        (altType ? word(name.substr(0,altPos), false) : name);

    const entry* eptr =
    (
        (dictptr != nullptr)
      ? dictptr->findScoped(lookupName, keyType::LITERAL_RECURSIVE)
      : nullptr
    );

    string value;
    if (eptr)
    {
        value = entryToString(eptr, allowSubDict);
    }
    else if (allowEnv || dictptr == nullptr)
    {
        value = Foam::getEnv(lookupName);
    }

    if (value.empty() ? (altType == '-') : (altType == '+'))
    {
        // Not found or empty:  use ":-" alternative value
        // Found and not empty: use ":+" alternative value
        value = name.substr(altPos + 2);
    }

    if (!allowEmpty && value.empty())
    {
        if (dictptr != nullptr)
        {
            auto& err =
                FatalIOErrorInFunction(*dictptr)
                    << "Cannot find dictionary entry ";

            if (allowEnv)
            {
                err << "or environment ";
            }

            err << "variable '" << lookupName << "'" << nl
                << exit(FatalIOError);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown variable '" << lookupName << "'" << nl
                << exit(FatalError);
        }
    }

    return value;
}


// Recursively expands (dictionary or environment) variable
// starting at index in string. Updates index.
//
// String:    "abc ${var} def",
// Receive:   "var} def"
//
// String:    "abc ${{expr}} def"
// Receive:   "{expr}} def"
//
// On return, the index will be adjust to be AFTER the closing '}'
static Foam::string recursiveExpand
(
    const std::string& s,
    std::string::size_type& index,
    const dictionary* dictptr,
    const bool allowEnv,
    const bool allowEmpty,
    const bool allowSubDict
)
{
    ///Info<< "process:" << index << "=" << s.substr(index) << endl;

    // Track ${{ expr }} expressions
    const bool isExpr = (index < s.size() && s[index] == '{');

    if (isExpr)
    {
        ++index;
    }

    // Initially called for a ${variable}, not ${{expr}}
    bool isVar = !isExpr;

    string out;

    for (/*nil*/; index < s.size(); ++index)
    {
        ///Info<< "remaining:" << index << "=" << s.substr(index) << endl;
        if (s[index] == '$')
        {
            if (s[index+1] == '{')
            {
                // Recurse to parse variable name
                index += 2;

                string val =
                    recursiveExpand
                    (
                        s,
                        index,
                        dictptr,
                        allowEnv,
                        allowEmpty,
                        allowSubDict
                    );

                out.append(val);    // Append content

                ///Info<< "got:" << val << nl << "now:" << out << endl;

                // Already skipped past '}' terminator?
                if (s[index-1] == '}')
                {
                    --index;
                }
            }
            else if (validVariableChar(s[index+1]))
            {
                // A regular $var expansion without a surrounding {}.

                const auto varLen = findVariableLen(s, index);
                const word varName(s.substr(index+1, varLen), false);
                index += varLen;

                string val =
                    getVariable
                    (
                        varName,
                        dictptr,
                        allowEnv,
                        allowEmpty,
                        allowSubDict
                    );

                out.append(val);    // Append content
            }
            else
            {
                // Something like '$()', '$[]', etc - pass through
                out += s[index];    // Append char
            }
        }
        else if (s[index] == '}')
        {
            // Closing an expression or variable

            if (isExpr)
            {
                // Closes with '}}'
                ++index;        // Index past closing '}'

                if (s[index] == '}')
                {
                    ++index;    // Index past closing '}'
                }
                else if (dictptr != nullptr)
                {
                    // Missing '}'? - Warn/error/ignore
                    FatalIOErrorInFunction(*dictptr)
                        << "Expansion ${{ is missing a closing '}}'\n"
                        << exit(FatalIOError);
                }
                else
                {
                    FatalErrorInFunction
                        << "Expansion ${{ is missing a closing '}}'\n"
                        << exit(FatalError);
                }

                ///Info<< "eval <" << out << ">" << endl;

                // Even with allow empty, expressions may need content

                string val(stringOps::evaluate(out));
                stringOps::inplaceTrim(val);

                return val;
            }
            else if (isVar)
            {
                // Variable - closes with '}'

                ++index;  // Index past closing '}'

                return
                    getVariable
                    (
                        out,
                        dictptr,
                        allowEnv,
                        allowEmpty,
                        allowSubDict
                    );
            }
            else
            {
                // Stray '}'? - Leave on output

                out += s[index];  // append char
            }
        }
        else
        {
            out += s[index];   // append char
        }
    }

    return out;
}


static void expandString
(
    std::string& s,
    const dictionary* dictptr,
    const bool allowEnv,
    const bool allowEmpty,
    const bool allowSubDict,
    const char sigil
)
{
    std::string::size_type varBeg = 0;

    // Expand $VAR, ${VAR} or ${{EXPR}}
    // Repeat until nothing more is found
    while
    (
        (varBeg = s.find(sigil, varBeg)) != std::string::npos
     && varBeg < s.size()-1
    )
    {
        if (varBeg == 0 || s[varBeg-1] != '\\')
        {
            if (s[varBeg+1] == '{')
            {
                // Recursive variable expansion mode: '${' or '${{'
                const auto replaceBeg = varBeg;

                varBeg += 2;
                string varValue
                (
                    recursiveExpand
                    (
                        s,
                        varBeg,
                        dictptr,
                        allowEnv,
                        allowEmpty,
                        allowSubDict
                    )
                );

                s.replace(replaceBeg, varBeg - replaceBeg, varValue);
                varBeg = replaceBeg+varValue.size();
            }
            else if (validVariableChar(s[varBeg+1]))
            {
                // A regular $var expansion without surrounding {}.
                const auto varLen(findVariableLen(s, varBeg, sigil));
                const word varName(s.substr(varBeg+1, varLen), false);

                string varValue
                (
                    getVariable
                    (
                        varName,
                        dictptr,
                        allowEnv,
                        allowEmpty,
                        allowSubDict
                    )
                );

                s.replace(varBeg, varName.size()+1, varValue);
                varBeg += varValue.size();
            }
            else
            {
                ++varBeg;
            }
        }
        else
        {
            ++varBeg;
        }
    }

    expandLeading(s);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

std::string::size_type Foam::stringOps::count(const char* s, const char c)
{
    // or: size_t len = Foam::string::length(s);
    size_t len = (s ? strlen(s) : 0);
    return len ? std::count(s, (s + len), c) : 0;
}


std::string::size_type Foam::stringOps::count
(
    const std::string& s,
    const char c
)
{
    return std::count(s.cbegin(), s.cend(), c);
}


Foam::string Foam::stringOps::expand
(
    const std::string& s,
    const HashTable<string>& mapping,
    const char sigil
)
{
    string out(s);
    inplaceExpand(out, mapping);
    return out;
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const HashTable<string>& mapping,
    const char sigil
)
{
    std::string::size_type varBeg = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (varBeg = s.find(sigil, varBeg)) != std::string::npos
     && varBeg < s.size()-1
    )
    {
        if (varBeg == 0 || s[varBeg-1] != '\\')
        {
            // Find end of first occurrence
            std::string::size_type varEnd = varBeg;
            std::string::size_type delim = 0;

            // The type/position of the ":-" or ":+" alternative values
            int altType = 0;
            auto altPos = std::string::npos;

            if (s[varBeg+1] == '{')
            {
                varEnd = s.find('}', varBeg);
                delim = 1;

                // Check for ${parameter:-word} or ${parameter:+word}
                if (varEnd != std::string::npos)
                {
                    altPos = varBeg;
                    altType = findParameterAlternative(s, altPos, varEnd);
                }
            }
            else
            {
                varEnd += findVariableLen(s, varBeg, sigil);
            }

            if (varEnd == std::string::npos)
            {
                // Likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (varEnd == varBeg)
            {
                // Something like '$()', '$[]', etc - pass through
                ++varBeg;
            }
            else
            {
                const word varName
                (
                    s.substr
                    (
                        varBeg + 1 + delim,
                        (
                            (altPos == std::string::npos ? varEnd : altPos)
                          - varBeg - 2*delim
                        )
                    ),
                    false
                );

                std::string altValue;
                if (altPos != std::string::npos)
                {
                    // Had ":-" or ":+" alternative value
                    altValue = s.substr
                    (
                        altPos + 2,
                        varEnd - altPos - 2*delim
                    );
                }


                const auto fnd = mapping.cfind(varName);

                if (fnd.good() ? (altType == '+') : (altType == '-'))
                {
                    // Found and ":+" alternative
                    // Not-found and ":-" alternative

                    s.replace(varBeg, varEnd - varBeg + 1, altValue);
                    varBeg += altValue.size();
                }
                else if (fnd.good())
                {
                    // Found: use value
                    s.replace(varBeg, varEnd - varBeg + 1, *fnd);
                    varBeg += (*fnd).size();
                }
                else
                {
                    // Not-found: empty value
                    s.erase(varBeg, varEnd - varBeg + 1);
                }
            }
        }
        else
        {
            ++varBeg;
        }
    }
}


Foam::string Foam::stringOps::expand
(
    const std::string& s,
    const dictionary& dict,
    const char sigil
)
{
    string out(s);
    inplaceExpand(out, dict, sigil);
    return out;
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const dictionary& dict,
    const bool allowEnv,
    const bool allowEmpty,
    const bool allowSubDict,
    const char sigil
)
{
    expandString(s, &dict, allowEnv, allowEmpty, allowSubDict, sigil);
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const dictionary& dict,
    const char sigil
)
{
    // Allow everything, including subDict expansions
    // env=true, empty=true, subDict=true
    expandString(s, &dict, true, true, true, sigil);
}


Foam::string Foam::stringOps::expand
(
    const std::string& s,
    const bool allowEmpty
)
{
    string out(s);
    inplaceExpand(out, allowEmpty);
    return out;
}


void Foam::stringOps::inplaceExpand
(
    std::string& s,
    const bool allowEmpty
)
{
    // Expand without a dictionary context
    // allowEnv=true, allowSubDict=N/A
    expandString(s, nullptr, true, allowEmpty, false, '$');
}


bool Foam::stringOps::inplaceReplaceVar(std::string& s, const word& varName)
{
    if (s.empty() || varName.empty())
    {
        return false;
    }

    const string content(Foam::getEnv(varName));
    if (content.empty())
    {
        return false;
    }

    const auto i = s.find(content);
    if (i == std::string::npos)
    {
        return false;
    }

    s.replace(i, content.size(), string("${" + varName + "}"));
    return true;
}


Foam::string Foam::stringOps::trimLeft(const std::string& s)
{
    if (!s.empty())
    {
        std::string::size_type pos = 0;
        const auto end = s.length();

        while (pos < end && std::isspace(s[pos]))
        {
            ++pos;
        }

        if (pos)
        {
            return s.substr(pos);
        }
    }

    return s;
}


void Foam::stringOps::inplaceTrimLeft(std::string& s)
{
    if (!s.empty())
    {
        std::string::size_type pos = 0;
        const auto end = s.length();

        while (pos < end && std::isspace(s[pos]))
        {
            ++pos;
        }

        if (pos)
        {
            s.erase(0, pos);
        }
    }
}


Foam::string Foam::stringOps::trimRight(const std::string& s)
{
    if (!s.empty())
    {
        auto end = s.length();
        while (end && std::isspace(s[end-1]))
        {
            --end;
        }

        if (end < s.length())
        {
            return s.substr(0, end);
        }
    }

    return s;
}


void Foam::stringOps::inplaceTrimRight(std::string& s)
{
    if (!s.empty())
    {
        auto end = s.length();
        while (end && std::isspace(s[end-1]))
        {
            --end;
        }

        s.erase(end);
    }
}


std::pair<std::size_t, std::size_t>
Foam::stringOps::findTrim
(
    const std::string& s,
    std::size_t pos,
    std::size_t len
)
{
    size_t end = s.length();
    if (pos >= end)
    {
        pos = end;
    }
    else if (len != std::string::npos)
    {
        len += pos;

        if (len < end)
        {
            end = len;
        }
    }

    // Right = last
    while (pos < end && std::isspace(s[end-1]))
    {
        --end;
    }

    // Left = first
    while (pos < end && std::isspace(s[pos]))
    {
        ++pos;
    }

    return std::pair<std::size_t, std::size_t>(pos, end);
}


Foam::string Foam::stringOps::trim(const std::string& s)
{
    std::string::size_type pos = 0;
    std::string::size_type end = s.length();

    // Right
    while (pos < end && std::isspace(s[end-1]))
    {
        --end;
    }

    // Left
    while (pos < end && std::isspace(s[pos]))
    {
        ++pos;
    }

    return s.substr(pos, end-pos);
}


void Foam::stringOps::inplaceTrim(std::string& s)
{
    inplaceTrimRight(s);
    inplaceTrimLeft(s);
}


void Foam::stringOps::inplaceRemoveSpace(std::string& s)
{
    s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());
}


Foam::string Foam::stringOps::removeComments(const std::string& s)
{
    string out(s);
    inplaceRemoveComments(out);
    return out;
}


void Foam::stringOps::inplaceRemoveComments(std::string& s)
{
    const auto len = s.length();

    if (len < 2)
    {
        return;
    }

    std::string::size_type n = 0;

    for (std::string::size_type i = 0; i < len; ++i)
    {
        char c = s[i];

        if (n != i)
        {
            s[n] = c;
        }
        ++n;

        // The start of a C/C++ comment?
        if (c == '/')
        {
            ++i;

            if (i == len)
            {
                // No further characters
                break;
            }

            c = s[i];

            if (c == '/')
            {
                // C++ comment - remove to end-of-line

                --n;
                s[n] = '\n';

                // Backtrack to eliminate leading spaces,
                // up to the previous newline

                while (n && std::isspace(s[n-1]))
                {
                    --n;

                    if (s[n] == '\n')
                    {
                        break;
                    }

                    s[n] = '\n';
                }

                i = s.find('\n', ++i);

                if (i == std::string::npos)
                {
                    // Truncated - done
                    break;
                }

                ++n;  // Include newline in output
            }
            else if (c == '*')
            {
                // C comment - search for '*/'
                --n;
                i = s.find("*/", ++i, 2);

                if (i == std::string::npos)
                {
                    // Truncated - done
                    break;
                }

                ++i;  // Index past first of "*/", loop increment does the rest
            }
            else
            {
                // Not a C/C++ comment
                if (n != i)
                {
                    s[n] = c;
                }
                ++n;
            }
        }
    }

    s.erase(n);
}


Foam::string Foam::stringOps::lower(const std::string& s)
{
    string out;
    out.resize(s.length());

    std::transform(s.begin(), s.end(), out.begin(), ::tolower);
    return out;
}


void Foam::stringOps::inplaceLower(std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}


Foam::string Foam::stringOps::upper(const std::string& s)
{
    string out;
    out.resize(s.length());

    std::transform(s.begin(), s.end(), out.begin(), ::toupper);
    return out;
}


void Foam::stringOps::inplaceUpper(std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
}


void Foam::stringOps::writeWrapped
(
    OSstream& os,
    const std::string& str,
    const std::string::size_type width,
    const std::string::size_type indent,
    const bool escape
)
{
    // Disabled below some minimal lower limit
    if (width <= 8)
    {
        char c = 0;

        const auto len = str.size();

        for (std::string::size_type pos = 0; pos < len; ++pos)
        {
            c = str[pos];

            if (escape && c == '\\')
            {
                os << '\\';
            }
            os << c;
        }

        // Trailing newline for non-empty string and if still pending
        if (len && c != '\n')
        {
            os << '\n';
        }

        return;
    }


    // Normal case
    std::size_t pos = 0;
    const auto len = str.size();

    // Output leading newlines without any indention
    while (pos < len && str[pos] == '\n')
    {
        os << '\n';
        ++pos;
    }

    while (pos < len)
    {
        // Potential end point, break point and next point
        std::string::size_type endp  = pos + width;
        std::string::size_type breakp = str.find('\n', pos);
        std::string::size_type nextp = endp;

        if (std::string::npos != breakp && breakp < endp)
        {
            // Embedded line break
            endp = breakp;
            nextp = breakp + 1;  // Skip this newline in the next chunk

            // Trim trailing space
            while
            (
                (endp > pos)
             && (str[endp-1] == ' ' || str[endp-1] == '\t')
            )
            {
                --endp;
            }
        }
        else if (endp >= len)
        {
            // Can output the rest without any wrapping, no line-breaks
            nextp = endp = len;
        }
        else
        {
            // Find a good point to break the string
            // try to find space/tab, or use punctuation as a fallback

            breakp = nextp = endp;
            std::string::size_type punc = std::string::npos;

            // Backtrack to find whitespace
            bool foundBreak = false;
            while (breakp > pos)
            {
                --breakp;

                const char c = str[breakp];

                if (c == ' ' || c == '\t')
                {
                    foundBreak = true;
                    endp = breakp;
                    // Found a space, but continue loop anyhow
                    // (trims trailing space)
                }
                else if (foundBreak)
                {
                    // Non-whitespace encountered while consuming
                    // trailing space. We are done
                    break;
                }
                else
                {
                    // Potentially viable as last non-whitespace?
                    nextp = breakp;

                    // Remember if we see any punctuation characters
                    // - useful later as fallback
                    if (punc == std::string::npos)
                    {
                        switch (c)
                        {
                            // Break before the punctuation
                            case '(' : case '<' :
                            {
                                punc = breakp;
                                break;
                            }
                            // Break after the punctuation
                            case ')' : case '>' :
                            case ',' : case '.' :
                            case ':' : case ';' :
                            case '/' : case '|' :
                            {
                                punc = (breakp + 1);
                                break;
                            }
                        }
                    }
                }
            }

            if (!foundBreak)
            {
                // No whitespace breaks, but perhaps a punctuation break.
                // Otherwise can't do much else

                if (punc != std::string::npos)
                {
                    nextp = endp = punc;
                }
                else
                {
                    nextp = endp;
                }
            }
        }


        // Output
        // ~~~~~~
        // Indent subsequent lines.
        // - assuming the one was done prior to calling this routine.
        // - no extra indent if it will only have a newline

        if (pos && (pos < endp))
        {
            // Put indent
            for (std::string::size_type i = 0; i < indent; ++i)
            {
                os << ' ';
            }
        }

        while (pos < endp)
        {
            const char c = str[pos];

            if (escape && c == '\\')
            {
                os << '\\';
            }
            os << c;

            ++pos;
        }
        os << nl;

        pos = nextp;
    }
}


// ************************************************************************* //

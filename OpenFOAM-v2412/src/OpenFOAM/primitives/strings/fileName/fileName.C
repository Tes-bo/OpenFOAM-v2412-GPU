/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "fileName.H"
#include "wordRe.H"
#include "wordList.H"
#include "DynamicList.H"
#include "OSspecific.H"
#include "fileOperation.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::fileName::typeName = "fileName";

int Foam::fileName::debug(Foam::debug::debugSwitch(fileName::typeName, 0));

int Foam::fileName::allowSpaceInFileName
(
    #ifdef _WIN32
    1  // Windows: expect spaces to occur
    #else
    Foam::debug::infoSwitch("allowSpaceInFileName", 0)
    #endif
);

const Foam::fileName Foam::fileName::null;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// doClean:
//   - remove duplicate slashes, "/./" and "/../" components.
//
// checkValid:
//   - similar to stripInvalid (but silent)
//
// return True if the content changed
static bool cleanFileName
(
    std::string& str,
    const bool doClean,
    const bool checkValid
)
{
    const auto maxLen = str.length();
    std::string::size_type nChar = 0;

    // Preserve UNC \\server\path (windows)
    // - MS-windows only, but handle for other systems
    //   since there is no collision with this pattern
    if (maxLen > 2 && str[0] == '\\' && str[1] == '\\')
    {
        nChar += 2;
    }

    char prev = 0;
    auto top = std::string::npos;  // Not yet found
    bool changed = false;

    for (auto src = nChar; src < maxLen; /*nil*/)
    {
        char c = str[src++];

        // Treat raw backslash like a path separator.
        // There is no "normal" way for these to be there
        // (except for an OS that uses them), but can cause issues
        // when writing strings, shell commands etc.
        if (c == '\\')
        {
            c = '/';
            str[nChar] = c;
            changed = true;
        }
        else if (checkValid && !Foam::fileName::valid(c))
        {
            // Ignore invalid chars
            // Could explicitly allow space character or rely on
            // allowSpaceInFileName via fileName::valid()
            continue;
        }

        if (c == '/' && top == std::string::npos)
        {
            // Top-level slash not previously determined
            top = (src-1);
        }

        if (doClean && prev == '/')
        {
            // Repeated '/' - skip it
            if (c == '/')
            {
                continue;
            }

            // Could be "/./", "/../" or a trailing "/."
            if (c == '.')
            {
                // Trailing "/." - skip it
                if (src >= maxLen)
                {
                    break;
                }

                // Peek at the next character
                const char c1 = str[src];

                // Found "/./" - skip over it
                if (c1 == '/' || c1 == '\\')
                {
                    ++src;
                    continue;
                }

                // Trailing "/.." or intermediate "/../"
                if
                (
                    c1 == '.'
                 &&
                    (
                        src+1 >= maxLen
                     || str[src+1] == '/' || str[src+1] == '\\'
                    )
                )
                {
                    // Backtrack to find the parent directory
                    // Minimum of 3 characters:  '/x/../'
                    // Strip it, provided it is above the top point

                    std::string::size_type parent;
                    if
                    (
                        nChar > 2
                     && top != std::string::npos
                     && (parent = str.rfind('/', nChar-2)) != std::string::npos
                     && parent >= top
                    )
                    {
                        nChar = parent + 1;   // Retain '/' from the parent
                        src += 2;
                        continue;
                    }

                    // Bad resolution, eg 'abc/../../'
                    // Retain the sequence, but move the top to avoid it being
                    // considered a valid parent later
                    top = nChar + 2;
                }
            }
        }

        str[nChar++] = prev = c;
    }

    // Remove trailing '/'
    if (doClean && nChar > 1 && str[nChar-1] == '/')
    {
        --nChar;
    }

    str.erase(nChar);
    return changed || (nChar != maxLen);
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::fileName::clean(std::string& str)
{
    return cleanFileName(str, true, false);  // clean, checkValid = false
}


Foam::fileName Foam::fileName::validate
(
    const std::string& str,
    const bool doClean
)
{
    fileName out(str, false);           // copy, no stripping
    cleanFileName(out, doClean, true);  // checkValid = true
    return out;
}


Foam::fileName Foam::fileName::concat
(
    const std::string& s1,
    const std::string& s2,
    const char delim
)
{
    const auto n1 = s1.length();
    const auto n2 = s2.length();

    fileName out;
    out.reserve(n1 + n2 + 1);

    out += s1;

    if (n1 && n2 && s1.back() != delim && s2.front() != delim)
    {
        // Add delimiter
        out += delim;
    }

    out += s2;

    // Could also remove trailing '/', if desired.
    return out;
}


bool Foam::fileName::equals(const std::string& s1, const std::string& s2)
{
    // Do not use (s1 == s2) or s1.compare(s2) first since this would
    // potentially be doing the comparison twice.

    std::string::size_type i1 = 0;
    std::string::size_type i2 = 0;

    const auto n1 = s1.length();
    const auto n2 = s2.length();

    //Info<< "compare " << s1 << " == " << s2 << endl;
    while (i1 < n1 && i2 < n2)
    {
        //Info<< "check '" << s1[i1] << "' vs '" << s2[i2] << "'" << endl;

        if (s1[i1] != s2[i2])
        {
            return false;
        }

        // Increment to next positions and also skip repeated slashes
        do
        {
            ++i1;
        }
        while (s1[i1] == '/');

        do
        {
            ++i2;
        }
        while (s2[i2] == '/');
    }
    //Info<< "return: " << Switch(i1 == n1 && i2 == n2) << endl;

    // Equal if it made it all the way through both strings
    return (i1 == n1 && i2 == n2);
}


bool Foam::fileName::isBackup(const std::string& s)
{
    if (s.empty())
    {
        return false;
    }
    else if (s.back() == '~')
    {
        return true;
    }

    // Now check the extension
    auto dot = find_ext(s);

    if (dot == npos)
    {
        return false;
    }

    ++dot;

    return
    (
        !s.compare(dot, npos, "bak") || !s.compare(dot, npos, "BAK")
     || !s.compare(dot, npos, "old") || !s.compare(dot, npos, "save")
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileName::fileName(const UList<word>& list)
{
    size_type len = 0;
    for (const word& item : list)
    {
        len += 1 + item.length();   // Include space for '/' needed
    }
    reserve(len);

    for (const word& item : list)
    {
        if (item.length())
        {
            if (length()) operator+=('/');
            operator+=(item);
        }
    }
}


Foam::fileName::fileName(std::initializer_list<word> list)
{
    size_type len = 0;
    for (const word& item : list)
    {
        len += 1 + item.length();   // Include space for '/' needed
    }
    reserve(len);

    for (const word& item : list)
    {
        if (item.length())
        {
            if (length()) operator+=('/');
            operator+=(item);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName::Type Foam::fileName::type
(
    bool followLink,
    bool checkGzip
) const
{
    Type t = ::Foam::type(*this, followLink);

    if (checkGzip && (Type::UNDEFINED == t) && size())
    {
        // Also check for gzip file?
        t = ::Foam::type(*this + ".gz", followLink);
    }

    return t;
}


Foam::fileName& Foam::fileName::toAbsolute()
{
    fileName& f = *this;

    if (!f.isAbsolute())
    {
        f = Foam::cwd()/f;
    }

    f.clean();  // Remove unneeded ".."

    return *this;
}


bool Foam::fileName::clean()
{
    return fileName::clean(*this);
}


std::string Foam::fileName::stem(const std::string& str)
{
    auto beg = str.rfind('/');
    auto dot = str.rfind('.');

    if (beg == npos)
    {
        beg = 0;
    }
    else
    {
        ++beg;
    }

    if (dot != npos && dot <= beg)
    {
        dot = npos;
    }

    if (dot == npos)
    {
        return str.substr(beg);
    }

    return str.substr(beg, dot - beg);
}


Foam::fileName& Foam::fileName::replace_name(const word& newName)
{
    const auto len = newName.length();

    if (len)
    {
        auto beg = rfind('/');

        if (beg == npos)
        {
            fileName::assign(newName);
        }
        else
        {
            ++beg;
            replace(beg, length()-beg, newName);
        }
    }

    return *this;
}


Foam::fileName Foam::fileName::relative
(
    const fileName& parent,
    const bool caseTag
) const
{
    const auto top = parent.length();
    const fileName& f = *this;

    // Everything after "parent/xxx/yyy" -> "xxx/yyy"
    //
    // case-relative:
    //     "parent/xxx/yyy" -> "<case>/xxx/yyy"
    //     "parent/constant/xxx/yyy" -> "<constant>/xxx/yyy"
    //     "parent/system/xxx/yyy" -> "<system>/xxx/yyy"
    //
    // as per stringOps::inplaceExpand()

    if
    (
        top && (f.length() > (top+1)) && f[top] == '/'
     && f.starts_with(parent)
    )
    {
        if (caseTag)
        {
            const auto trailing = f.find('/', top+1);

            if (npos != trailing)
            {
                switch (trailing-top-1)
                {
                    case 6:  // "system"
                    {
                        if (!compare((top+1), 6, "system"))
                        {
                            return "<system>"/f.substr(trailing+1);
                        }
                        break;
                    }
                    case 8:  // "constant"
                    {
                        if (!compare((top+1), 8, "constant"))
                        {
                            return "<constant>"/f.substr(trailing+1);
                        }
                        break;
                    }
                }
            }
            return "<case>"/f.substr(top+1);
        }
        else
        {
            return f.substr(top+1);
        }
    }
    else if (caseTag && f.length() && !f.isAbsolute())
    {
        const auto trailing = f.find('/');

        if (npos != trailing)
        {
            switch (trailing)
            {
                case 6:  // "system"
                {
                    if (!compare(0, 6, "system"))
                    {
                        return "<system>"/f.substr(trailing+1);
                    }
                    break;
                }
                case 8:  // "constant"
                {
                    if (!compare(0, 8, "constant"))
                    {
                        return "<constant>"/f.substr(trailing+1);
                    }
                    break;
                }
            }
        }
        return "<case>"/f;
    }

    return f;
}


Foam::wordList Foam::fileName::components(const char delim) const
{
    const auto parsed = stringOps::split<string>(*this, delim);

    wordList words(parsed.size());

    label i = 0;
    for (const auto& sub : parsed)
    {
        // Could easily filter out '.' here too
        words[i] = sub.str();
        ++i;
    }

    // As a plain wordList
    return words;
}


Foam::word Foam::fileName::component
(
    const size_type cmpt,
    const char delim
) const
{
    const auto parsed = stringOps::split<string>(*this, delim);

    if (parsed.size())
    {
        if (cmpt == std::string::npos)
        {
            return parsed.back().str();
        }
        else if (cmpt < parsed.size())
        {
            return parsed[cmpt].str();
        }
    }

    return word();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::fileName& Foam::fileName::operator/=(const string& other)
{
    fileName& s = *this;

    if (s.size())
    {
        if (other.size())
        {
            // Two non-empty strings: can concatenate

            if (s.back() != '/' && other.front() != '/')
            {
                s += '/';
            }

            s += other;
        }
    }
    else if (other.size())
    {
        // The first string is empty
        s = other;
    }

    return *this;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::fileName Foam::operator/(const string& s1, const string& s2)
{
    if (s1.length())
    {
        if (s2.length())
        {
            // Two non-empty strings: can concatenate

            if (s1.back() == '/' || s2.front() == '/')
            {
                return fileName(s1 + s2);
            }
            else
            {
                return fileName(s1 + '/' + s2);
            }
        }

        // The second string was empty
        return s1;
    }

    if (s2.length())
    {
        // The first string is empty
        return s2;
    }

    // Both strings are empty
    return fileName();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::search(const word& file, const fileName& directory)
{
    // Search current directory for the file
    for
    (
        const fileName& item
      : fileHandler().readDir(directory, fileName::FILE)
    )
    {
        if (item == file)
        {
            return directory/item;
        }
    }

    // If not found search each of the sub-directories
    for
    (
        const fileName& item
      : fileHandler().readDir(directory, fileName::DIRECTORY)
    )
    {
        fileName path = search(file, directory/item);
        if (!path.empty())
        {
            return path;
        }
    }

    return fileName();
}


// ************************************************************************* //

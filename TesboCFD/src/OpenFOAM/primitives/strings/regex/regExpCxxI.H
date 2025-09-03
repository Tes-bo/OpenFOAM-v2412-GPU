/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline std::regex::flag_type Foam::regExpCxx::syntax()
{
    // 0 = extended, 1 = ECMAScript
    return
    (
        regExpCxx::grammar
      ? std::regex::ECMAScript
      : std::regex::extended
    );
}


inline bool Foam::regExpCxx::is_meta(const char c) noexcept
{
    return
    (
        (c == '.')                           // any character
     || (c == '*' || c == '+' || c == '?')   // quantifiers
     || (c == '(' || c == ')' || c == '|')   // grouping/branching
     || (c == '[' || c == ']')               // range
    );
}


inline bool Foam::regExpCxx::is_meta
(
    const std::string& str,
    const char quote
)
{
    bool escaped = false;
    for (const char c : str)
    {
        if (quote && c == quote)
        {
            escaped = !escaped;  // toggle state
        }
        else if (escaped)
        {
            escaped = false;
        }
        else if (is_meta(c))
        {
            return true;
        }
    }
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::regExpCxx::regExpCxx()
:
    re_(),
    ctrl_(0)
{}


inline Foam::regExpCxx::regExpCxx(const regExpCxx& rgx)
:
    re_(rgx.re_),
    ctrl_(rgx.ctrl_)
{}


inline Foam::regExpCxx::regExpCxx(regExpCxx&& rgx) noexcept
:
    re_(std::move(rgx.re_)),
    ctrl_(rgx.ctrl_)
{
    rgx.ctrl_ = 0;
}


inline Foam::regExpCxx::regExpCxx
(
    const char* pattern,
    const bool ignoreCase
)
:
    re_(),
    ctrl_(0)
{
    set(pattern, ignoreCase);
}


inline Foam::regExpCxx::regExpCxx
(
    const std::string& pattern,
    const bool ignoreCase
)
:
    re_(),
    ctrl_(0)
{
    set(pattern, ignoreCase);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::regExpCxx::empty() const noexcept
{
    return !ctrl_;
}


inline bool Foam::regExpCxx::exists() const noexcept
{
    return ctrl_;
}


inline bool Foam::regExpCxx::negated() const noexcept
{
    return (ctrl_ == ctrlType::NEGATED);
}


inline bool Foam::regExpCxx::negate(bool on) noexcept
{
    bool old(ctrl_ == ctrlType::NEGATED);

    if (on)
    {
        if (ctrl_)
        {
            ctrl_ = ctrlType::NEGATED;
        }
    }
    else if (old)
    {
        ctrl_ = ctrlType::NORMAL;
    }

    return old;
}


inline unsigned Foam::regExpCxx::ngroups() const
{
    // Groups only make sense for regular (not negated) matching
    return ctrl_ == ctrlType::NORMAL ? re_.mark_count() : 0;
}


inline bool Foam::regExpCxx::nocase() const
{
    return ctrl_ && ((re_.flags() & std::regex::icase) == std::regex::icase);
}


inline bool Foam::regExpCxx::clear()
{
    if (ctrl_)
    {
        re_.assign("");
        ctrl_ = 0;

        return true;
    }

    return false;
}


inline void Foam::regExpCxx::swap(regExpCxx& rgx)
{
    if (this != &rgx)
    {
        // Self-swap is a no-op
        re_.swap(rgx.re_);
        std::swap(ctrl_, rgx.ctrl_);
    }
}


inline bool Foam::regExpCxx::set(const char* pattern, bool ignoreCase)
{
    // Silently handle nullptr
    return set_pattern
    (
        pattern,
        (pattern ? std::char_traits<char>::length(pattern) : 0),
        ignoreCase
    );
}


inline bool Foam::regExpCxx::set(const std::string& pattern, bool ignoreCase)
{
    return set_pattern
    (
        pattern.data(),
        pattern.length(),
        ignoreCase
    );
}


inline std::string::size_type
Foam::regExpCxx::find(const std::string& text) const
{
    // Find with negated is probably not very reliable...
    if (!ctrl_)
    {
        // Undefined: never matches
        return std::string::npos;
    }
    else if (text.empty())
    {
        if (ctrl_ == ctrlType::NEGATED)
        {
            return 0;  // No match - pretend it starts at position 0
        }
        else
        {
            return std::string::npos;
        }
    }
    else
    {
        std::smatch mat;

        const bool ok = std::regex_search(text, mat, re_);

        if (ctrl_ == ctrlType::NEGATED)
        {
            if (!ok)
            {
                return 0;  // No match - claim that is starts at position 0
            }
        }
        else if (ok)
        {
            return mat.position(0);
        }
    }

    return std::string::npos;  // False
}


inline bool Foam::regExpCxx::search(const std::string& text) const
{
    if (!ctrl_)
    {
        // Undefined: never matches
        return false;
    }

    const bool ok = (!text.empty() && std::regex_search(text, re_));

    return (ctrl_ == ctrlType::NEGATED) ? !ok : ok;
}


inline bool Foam::regExpCxx::match(const std::string& text) const
{
    if (!ctrl_)
    {
        // Undefined: never matches
        return false;
    }

    const bool ok = (!text.empty() && std::regex_match(text, re_));

    return (ctrl_ == ctrlType::NEGATED) ? !ok : ok;
}


inline bool Foam::regExpCxx::match
(
    const std::string& text,
    std::smatch& matches
) const
{
    // Probably does not make sense for negated pattern...
    if (negated())
    {
        // clear: matches = std::smatch();
        return match(text);
    }

    return std::regex_match(text, matches, re_);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline bool Foam::regExpCxx::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::regExpCxx::operator=(const regExpCxx& rgx)
{
    if (this != &rgx)
    {
        // Self-assignment is a no-op
        re_ = rgx.re_;
        ctrl_ = rgx.ctrl_;
    }
}


inline void Foam::regExpCxx::operator=(regExpCxx&& rgx)
{
    if (this != &rgx)
    {
        // Self-assignment is a no-op
        clear();
        swap(rgx);
    }
}


inline void Foam::regExpCxx::operator=(const char* pattern)
{
    set(pattern);
}


inline void Foam::regExpCxx::operator=(const std::string& pattern)
{
    set(pattern);
}


// ************************************************************************* //

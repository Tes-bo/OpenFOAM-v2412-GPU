/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::wordRe::valid(const char c)
{
    // Also accept '{' and '}' (for regex grouping?)
    return (word::valid(c) || c == '{' || c == '}');
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wordRe::wordRe()
:
    word(),
    regexPtr_(nullptr)
{}


inline Foam::wordRe::wordRe(const wordRe& str)
:
    word(static_cast<const word&>(str))
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(wordRe&& str)
:
    word(std::move(static_cast<word&>(str))),
    regexPtr_(str.regexPtr_.release())
{}


inline Foam::wordRe::wordRe(const word& str)
:
    word(str),
    regexPtr_(nullptr)
{}


inline Foam::wordRe::wordRe(word&& str)
:
    word(std::move(str)),
    regexPtr_(nullptr)
{}


inline Foam::wordRe::wordRe(const std::string& str, const compOption opt)
:
    word(str, false),  // No stripping
    regexPtr_(nullptr)
{
    if (opt != wordRe::LITERAL)
    {
        compile(opt);
    }
}


inline Foam::wordRe::wordRe(const char* str, const compOption opt)
:
    word(str, false),  // No stripping
    regexPtr_(nullptr)
{
    if (opt != wordRe::LITERAL)
    {
        compile(opt);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRe::isLiteral() const noexcept
{
    return !bool(regexPtr_);
}


inline bool Foam::wordRe::isPattern() const noexcept
{
    return bool(regexPtr_);
}


inline bool Foam::wordRe::compile(const compOption opt)
{
    if (opt != wordRe::LITERAL)
    {
        bool comp = false;

        if (opt & wordRe::REGEX)
        {
            comp = true;
        }
        else if (opt & wordRe::DETECT)
        {
            comp = regExp::is_meta(*this) || !string::valid<word>(*this);
        }
        else if (opt & wordRe::ICASE)
        {
            comp = true;
        }

        if (comp)
        {
            if (!regexPtr_)
            {
                regexPtr_.reset(new Foam::regExp());
            }

            if (!regexPtr_->set(*this, (opt & wordRe::ICASE)))
            {
                // Compilation failed
                regexPtr_.reset(nullptr);
            }

            return bool(regexPtr_);
        }
    }

    // Fall-through behaviour - not a regex
    regexPtr_.reset(nullptr);
    return false;
}


inline bool Foam::wordRe::compile()
{
    if (!regexPtr_)
    {
        regexPtr_.reset(new Foam::regExp());
    }

    if (!regexPtr_->set(*this))
    {
        // Compilation failed
        regexPtr_.reset(nullptr);
    }

    return bool(regexPtr_);
}


inline void Foam::wordRe::uncompile()
{
    regexPtr_.reset(nullptr);
}


inline void Foam::wordRe::uncompile(bool adjust)
{
    // Only strip when debug is active (potentially costly operation)
    if (adjust && isPattern() && word::debug)
    {
        string::stripInvalid<word>(*this);
    }
    regexPtr_.reset(nullptr);
}


inline void Foam::wordRe::clear()
{
    word::clear();
    regexPtr_.reset(nullptr);
}


inline bool Foam::wordRe::match(const std::string& text, bool literal) const
{
    if (!literal && regexPtr_)
    {
        return regexPtr_->match(text);  // Match as regex
    }

    return !compare(text);  // Compare as literal
}


inline void Foam::wordRe::set(const std::string& str, const compOption opt)
{
    assign(str);
    compile(opt);
}


inline void Foam::wordRe::set(const char* str, const compOption opt)
{
    // No nullptr protection here
    assign(str);
    compile(opt);
}


inline void Foam::wordRe::swap(wordRe& str)
{
    if (this == &str)
    {
        return;  // Self-swap is a no-op
    }

    word::swap(static_cast<word&>(str));
    regexPtr_.swap(str.regexPtr_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wordRe::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::wordRe::operator=(const wordRe& str)
{
    if (this == &str)
    {
        return;  // Self-assignment is a no-op
    }

    assign(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        regexPtr_.reset(nullptr);
    }
}


inline void Foam::wordRe::operator=(const word& str)
{
    assign(str);
    regexPtr_.reset(nullptr);
}


inline void Foam::wordRe::operator=(const string& str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const std::string& str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const char* str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(wordRe&& str)
{
    if (this == &str)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    swap(str);
}


// ************************************************************************* //

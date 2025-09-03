/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::keyType::valid(const char c)
{
    // Also accept '{' and '}' (for regex grouping?)
    return (word::valid(c) || c == '{' || c == '}');
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::keyType::keyType()
:
    word(),
    type_(option::LITERAL)
{}


inline Foam::keyType::keyType(const word& str)
:
    word(str),
    type_(option::LITERAL)
{}


inline Foam::keyType::keyType(word&& str)
:
    word(std::move(str)),
    type_(option::LITERAL)
{}


inline Foam::keyType::keyType(const string& str)
:
    word(str, false),  // No stripping
    type_(option::REGEX)
{}


inline Foam::keyType::keyType(string&& str)
:
    word(std::move(str), false),  // No stripping
    type_(option::REGEX)
{}


inline Foam::keyType::keyType(const std::string& str, option opt)
:
    word(str, false),  // No stripping
    type_(option(opt & 0x0F))
{}


inline Foam::keyType::keyType(std::string&& str, option opt)
:
    word(std::move(str), false),  // No stripping
    type_(option(opt & 0x0F))
{}


inline Foam::keyType::keyType(const char* str, option opt)
:
    word(str, false),  // No stripping
    type_(option(opt & 0x0F))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::isLiteral() const noexcept
{
    return !(type_ & option::REGEX);
}


inline bool Foam::keyType::isPattern() const noexcept
{
    return (type_ & option::REGEX);
}


inline void Foam::keyType::setType(option opt, bool adjust)
{
    opt = option(opt & 0x0F);

    if (type_ != opt)
    {
        // Only strip when debug is active (potentially costly operation)
        if (adjust && isPattern() && word::debug)
        {
            string::stripInvalid<word>(*this);
        }

        type_ = opt;
    }
}


inline bool Foam::keyType::compile() noexcept
{
    type_ = option::REGEX;
    return true;
}


inline void Foam::keyType::uncompile() noexcept
{
    type_ = option::LITERAL;
}


inline void Foam::keyType::uncompile(bool adjust)
{
    // Only strip when debug is active (potentially costly operation)
    if (adjust && isPattern() && word::debug)
    {
        string::stripInvalid<word>(*this);
    }

    type_ = option::LITERAL;
}


inline void Foam::keyType::clear()
{
    word::clear();
    type_ = option::LITERAL;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::keyType::operator()(const std::string& text) const
{
    return match(text);  // Use smart match
}


inline void Foam::keyType::operator=(const keyType& str)
{
    if (this == &str)
    {
        return;  // Self-assignment is a no-op
    }

    assign(str); // Bypasses char checking
    type_ = str.type_;
}


inline void Foam::keyType::operator=(keyType&& str)
{
    if (this == &str)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    swap(str);
}


inline void Foam::keyType::operator=(const char* str)
{
    assign(str); // Bypasses char checking
    type_ = option::LITERAL;
}


inline void Foam::keyType::operator=(const word& str)
{
    assign(str); // Bypasses char checking
    type_ = option::LITERAL;
}


inline void Foam::keyType::operator=(const string& str)
{
    assign(str); // Bypasses char checking
    type_ = option::REGEX;
}


// ************************************************************************* //

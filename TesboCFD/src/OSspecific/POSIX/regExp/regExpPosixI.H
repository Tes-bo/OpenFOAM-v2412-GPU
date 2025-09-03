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

#include <utility>  // For std::swap

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::regExpPosix::regExpPosix() noexcept
:
    preg_(nullptr),
    ctrl_(0)
{}


inline Foam::regExpPosix::regExpPosix(regExpPosix&& rgx) noexcept
:
    preg_(rgx.preg_),
    ctrl_(rgx.ctrl_)
{
    rgx.preg_ = nullptr;
    rgx.ctrl_ = 0;
}


inline Foam::regExpPosix::regExpPosix
(
    const char* pattern,
    const bool ignoreCase
)
:
    preg_(nullptr),
    ctrl_(0)
{
    set(pattern, ignoreCase);
}


inline Foam::regExpPosix::regExpPosix
(
    const std::string& pattern,
    const bool ignoreCase
)
:
    preg_(nullptr),
    ctrl_(0)
{
    set(pattern, ignoreCase);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::regExpPosix::~regExpPosix()
{
    clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::regExpPosix::empty() const noexcept
{
    return !preg_;
}


inline bool Foam::regExpPosix::exists() const noexcept
{
    return preg_;
}


inline bool Foam::regExpPosix::negated() const noexcept
{
    return (ctrl_ == ctrlType::NEGATED);
}


inline bool Foam::regExpPosix::negate(bool on) noexcept
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


inline unsigned Foam::regExpPosix::ngroups() const
{
    return (preg_ && ctrl_ == ctrlType::NORMAL) ? preg_->re_nsub : 0;
}


inline bool Foam::regExpPosix::set(const char* pattern, bool ignoreCase)
{
    // Silently handle nullptr
    return set_pattern
    (
        pattern,
        (pattern ? std::char_traits<char>::length(pattern) : 0),
        ignoreCase
    );
}


inline bool Foam::regExpPosix::set(const std::string& pattern, bool ignoreCase)
{
    return set_pattern
    (
        pattern.data(),
        pattern.length(),
        ignoreCase
    );
}


inline bool Foam::regExpPosix::search(const std::string& text) const
{
    if (!ctrl_)
    {
        return false;
    }
    else if (text.empty())
    {
        return (ctrl_ == ctrlType::NEGATED);
    }

    return std::string::npos != find(text);
}


inline void Foam::regExpPosix::swap(regExpPosix& rgx)
{
    if (this != &rgx)
    {
        // Self-swap is a no-op
        std::swap(preg_, rgx.preg_);
        std::swap(ctrl_, rgx.ctrl_);
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline bool Foam::regExpPosix::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::regExpPosix::operator=(regExpPosix&& rgx)
{
    if (this != &rgx)
    {
        // Self-assignment is a no-op
        clear();
        swap(rgx);
    }
}


inline void Foam::regExpPosix::operator=(const char* pattern)
{
    set(pattern);
}


inline void Foam::regExpPosix::operator=(const std::string& pattern)
{
    set(pattern);
}


// ************************************************************************* //

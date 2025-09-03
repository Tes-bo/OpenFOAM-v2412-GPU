/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline Foam::label Foam::wordRes::first_match
(
    const UList<wordRe>& selectors,
    const std::string& text,
    const bool literal
)
{
    label index = 0;
    for (const wordRe& select : selectors)
    {
        if (select.match(text, literal))
        {
            return index;
        }
        ++index;
    }

    return -1;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline bool Foam::wordRes::match
(
    const UList<wordRe>& selectors,
    const std::string& text,
    bool literal
)
{
    for (const wordRe& select : selectors)
    {
        if (select.match(text, literal))
        {
            return true;
        }
    }
    return false;
}


inline Foam::wordRe::compOption Foam::wordRes::matched
(
    const UList<wordRe>& selectors,
    const std::string& text
)
{
    auto retval(wordRe::compOption::UNKNOWN);

    for (const wordRe& select : selectors)
    {
        if (select.isLiteral())
        {
            if (select.match(text, true))
            {
                return wordRe::compOption::LITERAL;
            }
        }
        else if
        (
            // Only match regex once
            retval == wordRe::compOption::UNKNOWN
         && select.match(text, false)
        )
        {
            retval = wordRe::compOption::REGEX;
        }
    }

    return retval;
}


template<class StringType>
inline Foam::labelList Foam::wordRes::matching
(
    const wordRe& select,
    const UList<StringType>& input,
    const bool invert
)
{
    if (select.empty() && !invert)
    {
        return labelList();
    }

    const label len = input.size();

    labelList indices(len);

    label count = 0;
    for (label i = 0; i < len; ++i)
    {
        if (select.match(input[i]) ? !invert : invert)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


template<class StringType>
inline Foam::labelList Foam::wordRes::matching
(
    const UList<wordRe>& selectors,
    const UList<StringType>& input,
    const bool invert
)
{
    if (selectors.empty() && !invert)
    {
        return labelList();
    }

    const label len = input.size();

    labelList indices(len);

    label count = 0;
    for (label i = 0; i < len; ++i)
    {
        if (wordRes::match(selectors, input[i]) ? !invert : invert)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRes::match(const std::string& text, bool literal) const
{
    return wordRes::match(*this, text, literal);
}


inline Foam::wordRe::compOption
Foam::wordRes::matched(const std::string& text) const
{
    return wordRes::matched(*this, text);
}


template<class StringType>
inline Foam::labelList Foam::wordRes::matching
(
    const UList<StringType>& input,
    const bool invert
) const
{
    return wordRes::matching(*this, input, invert);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wordRes::operator()(const std::string& text) const
{
    return wordRes::match(*this, text);
}


// * * * * * * * * * * * * * * * * Functors  * * * * * * * * * * * * * * * * //

inline Foam::wordRes::matcher::matcher
(
    const UList<wordRe>& selectors
) noexcept
:
    select_(selectors)
{}


inline Foam::wordRes::filter::filter
(
    const UList<wordRe>& allow,
    const UList<wordRe>& deny
) noexcept
:
    allow_(allow),
    deny_(deny)
{}


inline bool Foam::wordRes::matcher::empty() const noexcept
{
    return select_.empty();
}

inline bool Foam::wordRes::filter::empty() const noexcept
{
    return (allow_.empty() && deny_.empty());
}


inline bool Foam::wordRes::matcher::operator()(const std::string& text) const
{
    return wordRes::match(select_, text);
}


inline bool Foam::wordRes::filter::operator()(const std::string& text) const
{
    if (allow_.empty())
    {
        // No allow specified, so accept everything that is NOT blocked
        return (deny_.empty() || !wordRes::match(deny_, text));
    }
    else if (deny_.empty())
    {
        // Nothing blocked, apply accept filter
        return wordRes::match(allow_, text);
    }
    else
    {
        // Both accept and deny filters, need to search more carefully
        const auto result = wordRes::matched(allow_, text);

        return
        (
            result == wordRe::LITERAL
          ? true
          :
            (
                result == wordRe::REGEX
             && !wordRes::match(deny_, text)
            )
        );
    }
}


// ************************************************************************* //

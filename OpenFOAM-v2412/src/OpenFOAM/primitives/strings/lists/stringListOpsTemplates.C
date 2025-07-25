/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class UnaryMatchPredicate, class StringType>
Foam::label Foam::firstMatchingString
(
    const UnaryMatchPredicate& matcher,
    const UList<StringType>& input,
    const bool invert
)
{
    const label len = input.size();

    for (label i=0; i < len; ++i)
    {
        if (matcher(input[i]) ? !invert : invert)
        {
            return i;
        }
    }

    return -1;
}


template<class UnaryMatchPredicate, class StringType>
Foam::labelList Foam::findMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    const UList<StringType>& input,
    const bool invert
)
{
    const label len = input.size();

    labelList indices(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (matcher(input[i]) ? !invert : invert)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


template<class UnaryMatchPredicate, class StringListType>
StringListType Foam::subsetMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    const StringListType& input,
    const bool invert
)
{
    const label len = input.size();

    StringListType output(len);
    output.resize(len);   // Consistent sizing (eg, DynamicList)

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (matcher(input[i]) ? !invert : invert)
        {
            output[count] = input[i];
            ++count;
        }
    }
    output.resize(count);

    return output;
}


template<class UnaryMatchPredicate, class StringListType>
void Foam::inplaceSubsetMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    StringListType& input,
    const bool invert
)
{
    const label len = input.size();

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (matcher(input[i]) ? !invert : invert)
        {
            if (count != i)
            {
                input[count] = std::move(input[i]);
            }
            ++count;
        }
    }
    input.resize(count);
}


template<class StringListType, class AccessOp>
Foam::labelList Foam::stringListOps::findMatching
(
    const StringListType& input,
    const wordRes::filter& pred,
    AccessOp aop
)
{
    const label len = input.size();

    if (pred.empty())
    {
        // Accept all
        return identity(len);
    }

    labelList indices(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        const std::string& text = aop(input[i]);

        if (pred(text))
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


template<class StringListType, class AccessOp>
Foam::labelList Foam::stringListOps::findMatching
(
    const StringListType& input,
    const wordRes& allow,
    const wordRes& deny,
    AccessOp aop
)
{
    if (allow.empty() && deny.empty())
    {
        // Accept all
        return identity(input.size());
    }

    // Use combined accept/reject filter
    const wordRes::filter pred(allow, deny);

    return stringListOps::findMatching(input, pred, aop);
}


// ************************************************************************* //

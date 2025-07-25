/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class IntListType>
IntListType Foam::renumber
(
    const labelUList& oldToNew,
    const IntListType& input
)
{
    const label len = input.size();

    IntListType output(len);
    output.resize(len);     // Consistent sizing (eg, DynamicList)

    for (label i=0; i < len; ++i)
    {
        if (input[i] >= 0)
        {
            output[i] = oldToNew[input[i]];
        }
    }

    return output;
}


template<class IntListType>
void Foam::inplaceRenumber
(
    const labelUList& oldToNew,
    IntListType& input
)
{
    const label len = input.size();

    for (label i=0; i < len; ++i)
    {
        if (input[i] >= 0)
        {
            input[i] = oldToNew[input[i]];
        }
    }
}


template<class IntListType>
void Foam::inplaceRenumber
(
    const Map<label>& oldToNew,
    IntListType& input
)
{
    const label len = input.size();

    for (label i = 0; i < len; ++i)
    {
        const auto fnd = oldToNew.cfind(input[i]);

        if (fnd.good())
        {
            input[i] = fnd.val();
        }
    }
}


template<class ListType>
ListType Foam::reorder
(
    const labelUList& oldToNew,
    const ListType& input,
    const bool prune
)
{
    const label len = input.size();

    ListType output(len);
    output.resize(len);     // Consistent sizing (eg, DynamicList)

    label maxIdx = -1;      // For pruning: The new size = maxIdx+1
    for (label i=0; i < len; ++i)
    {
        const label newIdx = oldToNew[i];
        if (newIdx >= 0)
        {
            // Could enforce (newIdx < len)
            // ... or just rely on FULLDEBUG from UList

            output[newIdx] = input[i];

            if (maxIdx < newIdx)
            {
                maxIdx = newIdx;
            }
        }
        else if (!prune)
        {
            output[i] = input[i];
        }
    }

    if (prune)
    {
        output.resize(maxIdx+1);
    }

    return output;
}


template<class ListType>
void Foam::inplaceReorder
(
    const labelUList& oldToNew,
    ListType& inputOutput,
    const bool prune
)
{
    // NOTE: cannot use std::move() since we have no guarantee that
    // the oldToNew map is unique (ie, shuffle)

    // Use const reference to ensure we obtain the proper operator[]
    // on lazy lists (eg, List<bool>, PackedList)

    const ListType& input = inputOutput;
    const label len = input.size();

    ListType output(len);
    output.resize(len);     // Consistent sizing (eg, DynamicList)

    label maxIdx = -1;      // For pruning: The new size = maxIdx+1
    for (label i=0; i < len; ++i)
    {
        const label newIdx = oldToNew[i];
        if (newIdx >= 0)
        {
            // Could enforce (newIdx < len)
            // ... or just rely on FULLDEBUG from UList

            output[newIdx] = input[i];

            if (maxIdx < newIdx)
            {
                maxIdx = newIdx;
            }
        }
        else if (!prune)
        {
            output[i] = input[i];
        }
    }

    if (prune)
    {
        output.resize(maxIdx+1);
    }

    inputOutput.transfer(output);
}


template<unsigned Width>
Foam::PackedList<Width> Foam::reorder
(
    const labelUList& oldToNew,
    const PackedList<Width>& input,
    const bool prune
)
{
    const label len = input.size();

    PackedList<Width> output(len);

    label maxIdx = -1;      // For pruning: The new size = maxIdx+1
    for (label i=0; i < len; ++i)
    {
        const auto& val = input.get(i);

        const label newIdx = oldToNew[i];

        if (newIdx >= 0)
        {
            // Could enforce (newIdx < len)
            // ... or just rely on FULLDEBUG from UList

            output.set(newIdx, val);

            if (maxIdx < newIdx)
            {
                maxIdx = newIdx;
            }
        }
        else if (!prune)
        {
            output.set(i, val);
        }
    }

    if (prune)
    {
        output.resize(maxIdx+1);
    }

    // Verify addresses (for movable refs)
    // Info<< "reordered in " << name(input.cdata()) << nl
    //     << "reordered out " << name(output.cdata()) << nl;

    return output;
}


template<unsigned Width>
void Foam::inplaceReorder
(
    const labelUList& oldToNew,
    PackedList<Width>& input,
    const bool prune
)
{
    input = Foam::reorder(oldToNew, input, prune);

    // Verify address (for movable refs)
    // Info<< "now have " << name(input.cdata()) << nl;
}


template<class Container>
void Foam::inplaceMapKey
(
    const labelUList& oldToNew,
    Container& input
)
{
    Container output(input.capacity());

    for (auto iter = input.begin(); iter != input.end(); ++iter)
    {
        const label oldIdx = iter.key();
        if (oldIdx >= 0)
        {
            // Could enforce (oldIdx < oldToNew.size())
            // ... or just rely on FULLDEBUG from UList

            output.insert(oldToNew[oldIdx], iter.val());
        }
    }

    input.transfer(output);
}


template<class Container>
Foam::label Foam::inplaceMapValue
(
    const labelUList& oldToNew,
    Container& input
)
{
    label nChanged = 0;

    for (auto iter = input.begin(); iter != input.end(); ++iter)
    {
        const label oldIdx = iter.val();
        if (oldIdx >= 0)
        {
            // Could enforce (oldIdx < oldToNew.size())
            // ... or just rely on FULLDEBUG from UList

            const label newIdx = oldToNew[oldIdx];

            if (oldIdx != newIdx)
            {
                iter.val() = newIdx;
                ++nChanged;
            }
        }
    }

    return nChanged;
}


template<class Container>
Foam::label Foam::inplaceMapValue
(
    const Map<label>& mapper,
    Container& input
)
{
    if (mapper.empty())
    {
        return 0;
    }

    label nChanged = 0;

    for (auto iter = input.begin(); iter != input.end(); ++iter)
    {
        label& value = iter.val();

        auto mapIter = mapper.find(value);
        if (mapIter.good() && value != mapIter.val())
        {
            value = mapIter.val();
            ++nChanged;
        }
    }

    return nChanged;
}


template<class T>
Foam::labelList Foam::duplicateOrder
(
    const UList<T>& input
)
{
    labelList order;
    duplicateOrder(input, order, typename UList<T>::less(input));
    return order;
}


template<class T>
void Foam::duplicateOrder
(
    const UList<T>& input,
    labelList& order
)
{
    duplicateOrder(input, order, typename UList<T>::less(input));
}


template<class T, class ListComparePredicate>
void Foam::duplicateOrder
(
    const UList<T>& input,
    labelList& order,
    const ListComparePredicate& comp
)
{
    if (input.size() < 2)
    {
        order.clear();
        return;
    }

    Foam::sortedOrder(input, order, comp);

    const label last = (order.size()-1);
    label count = 0;
    for (label i = 0; i < last; ++i)
    {
        if (input[order[i]] == input[order[i+1]])
        {
            order[count] = order[i];
            ++count;
        }
    }
    order.resize(count);
}


template<class T>
Foam::labelList Foam::uniqueOrder
(
    const UList<T>& input
)
{
    labelList order;
    uniqueOrder(input, order, typename UList<T>::less(input));
    return order;
}


template<class T>
void Foam::uniqueOrder
(
    const UList<T>& input,
    labelList& order
)
{
    uniqueOrder(input, order, typename UList<T>::less(input));
}


template<class T, class ListComparePredicate>
void Foam::uniqueOrder
(
    const UList<T>& input,
    labelList& order,
    const ListComparePredicate& comp
)
{
    Foam::sortedOrder(input, order, comp);

    if (order.size() > 1)
    {
        const label last = (order.size()-1);
        label count = 0;
        for (label i = 0; i < last; ++i)
        {
            if (input[order[i]] != input[order[i+1]])
            {
                order[count++] = order[i];
            }
        }
        order[count++] = order[last];
        order.resize(count);
    }
}


template<class T>
Foam::List<T> Foam::uniqueSort(const UList<T>& input)
{
    List<T> output(input);

    const label len = output.size();

    if (len > 1)
    {
        Foam::stableSort(output);

        label count = 0;
        for (label i = 1; i < len; ++i)
        {
            if (output[count] != output[i])
            {
                output[++count] = output[i];
            }
        }

        output.resize(count+1);
    }

    return output;
}


template<class ListType>
void Foam::inplaceUniqueSort(ListType& input)
{
    inplaceUniqueSort
    (
        input,
        typename UList<typename ListType::value_type>::less(input)
    );
}


template<class ListType, class ListComparePredicate>
void Foam::inplaceUniqueSort
(
    ListType& input,
    const ListComparePredicate& comp
)
{
    labelList order;
    uniqueOrder(input, order, comp);

    const label len = order.size();

    ListType output(len);
    output.resize(len);     // Consistent sizing (eg, DynamicList)

    for (label i=0; i < len; ++i)
    {
        output[i] = std::move(input[order[i]]);
    }

    input.transfer(output);
}


template<class BoolListType, class T>
Foam::List<T> Foam::subset
(
    const BoolListType& select,
    const UList<T>& input,
    const bool invert
)
{
    // Note: select can have a different size (eg, labelHashSet)

    const label len = input.size();

    List<T> output(len);

    label count = 0;

    for (label i=0; i < len; ++i)
    {
        if (select[i] ? !invert : invert)
        {
            output[count] = input[i];
            ++count;
        }
    }

    output.resize(count);

    return output;
}


template<class T>
Foam::List<T> Foam::subset
(
    const bitSet& select,
    const UList<T>& input,
    const bool invert
)
{
    const label len = input.size();

    List<T> output;

    label count = 0;

    if (!invert)
    {
        output.resize(select.count());

        for (const label i : select)
        {
            if (i >= len) break; // Avoid out of bounds (when select is longer)

            output[count] = input[i];
            ++count;
        }
    }
    else
    {
        const label outlen = (select.size() - select.count());
        output.resize(outlen);

        for (label i=0; i < len; ++i)
        {
            if (!select[i])
            {
                output[count] = input[i];
                ++count;
                if (count >= outlen) break;  // terminate early
            }
        }
    }

    output.resize(count);

    return output;
}


template<class BoolListType, class ListType>
void Foam::inplaceSubset
(
    const BoolListType& select,
    ListType& input,
    const bool invert
)
{
    // Note: select can have a different size (eg, labelHashSet)

    const label len = input.size();

    label count = 0;

    for (label i=0; i < len; ++i)
    {
        if (select[i] ? !invert : invert)
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


template<class ListType>
void Foam::inplaceSubset
(
    const bitSet& select,
    ListType& input,
    const bool invert
)
{
    label count = 0;

    if (!invert)
    {
        // Normal selection

        const label len = input.size();

        for (const label i : select)
        {
            if (i >= len) break;

            if (count != i)
            {
                input[count] = std::move(input[i]);
            }
            ++count;
        }
    }
    else
    {
        // Inverted selection

        const label outlen = (select.size() - select.count());

        const label len = min(input.size(), select.size());

        for (label i=0; i < len; ++i)
        {
            if (!select[i])
            {
                if (count != i)
                {
                    input[count] = std::move(input[i]);
                }
                ++count;
                if (count >= outlen) break;  // terminate early
            }
        }
    }

    input.resize(count);
}


template<class T, class UnaryPredicate>
Foam::List<T> Foam::subsetList
(
    const UList<T>& input,
    const UnaryPredicate& pred,
    const bool invert
)
{
    const label len = input.size();

    List<T> output(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (pred(input[i]) ? !invert : invert)
        {
            output[count] = input[i];
            ++count;
        }
    }

    output.resize(count);

    return output;
}


template<class ListType, class UnaryPredicate>
void Foam::inplaceSubsetList
(
    ListType& input,
    const UnaryPredicate& pred,
    const bool invert
)
{
    const label len = input.size();

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (pred(input[i]) ? !invert : invert)
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


template<class InputIntListType, class OutputIntListType>
void Foam::invertManyToMany
(
    const label len,
    const UList<InputIntListType>& input,
    List<OutputIntListType>& output
)
{
    // The output list sizes
    labelList sizes(len, Foam::zero{});

    for (const InputIntListType& sublist : input)
    {
        forAll(sublist, idx)
        {
            sizes[sublist[idx]]++;
        }
    }

    // Size output
    output.resize(len);
    forAll(sizes, outi)
    {
        output[outi].resize(sizes[outi]);
    }

    // Fill output
    sizes = 0;
    forAll(input, listi)
    {
        const InputIntListType& sublist = input[listi];

        forAll(sublist, idx)
        {
            const label outi = sublist[idx];

            output[outi][sizes[outi]++] = listi;
        }
    }
}


template<class ListType>
Foam::labelList Foam::findIndices
(
    const ListType& input,
    typename ListType::const_reference val,
    label start
)
{
    const label len = input.size();

    // Pass 1: count occurrences
    label count = 0;

    if (start >= 0)
    {
        for (label i = start; i < len; ++i)
        {
            if (input[i] == val)
            {
                if (!count) start = i;  // adjust start for second pass
                ++count;
            }
        }
    }

    labelList indices(count);

    // Pass 2: fill content
    if (count)
    {
        const label total = count;
        count = 0;
        for (label i = start; i < len; ++i)
        {
            if (input[i] == val)
            {
                indices[count] = i;
                if (++count == total)  // early termination
                {
                    break;
                }
            }
        }
    }

    return indices;
}


template<class ListType>
Foam::label Foam::findMin
(
    const ListType& input,
    label start
)
{
    const label len = input.size();

    if (start < 0 || start >= len)
    {
        return -1;
    }

    for (label i = start+1; i < len; ++i)
    {
        if (input[i] < input[start])
        {
            start = i;
        }
    }

    return start;
}


template<class ListType>
Foam::label Foam::findMax
(
    const ListType& input,
    label start
)
{
    const label len = input.size();

    if (start < 0 || start >= len)
    {
        return -1;
    }

    for (label i = start+1; i < len; ++i)
    {
        if (input[start] < input[i])
        {
            start = i;
        }
    }

    return start;
}


template<class ListType>
Foam::labelPair Foam::findMinMax
(
    const ListType& input,
    label start
)
{
    const label len = input.size();

    if (start < 0 || start >= len)
    {
        return labelPair(-1,-1);
    }

    label minIdx = start;
    label maxIdx = start;

    for (label i = start+1; i < len; ++i)
    {
        if (input[i] < input[minIdx])
        {
            minIdx = i;
        }
        if (input[maxIdx] < input[i])
        {
            maxIdx = i;
        }
    }

    return labelPair(minIdx, maxIdx);
}


template<class ListType>
Foam::label Foam::findSortedIndex
(
    const ListType& input,
    typename ListType::const_reference val,
    const label start
)
{
    label low = start;
    label high = input.size() - 1;

    if (start < 0 || start >= input.size())
    {
        return -1;
    }

    while (low <= high)
    {
        const label mid = (low + high)/2;

        if (val < input[mid])
        {
            high = mid - 1;
        }
        else if (input[mid] < val)
        {
            low = mid + 1;
        }
        else
        {
            return mid;
        }
    }

    return -1;
}


template<class ListType, class T, class ComparePredicate>
Foam::label Foam::findLower
(
    const ListType& input,
    const T& val,
    const label start,
    const ComparePredicate& comp
)
{
    label low = start;
    label high = input.size() - 1;

    if (start < 0 || start >= input.size())
    {
        return -1;
    }

    while ((high - low) > 1)
    {
        const label mid = (low + high)/2;

        if (comp(input[mid], val))
        {
            low = mid;
        }
        else
        {
            high = mid;
        }
    }

    if (comp(input[high], val))
    {
        return high;
    }
    else if (comp(input[low], val))
    {
        return low;
    }
    else
    {
        return -1;
    }
}


template<class ListType, class T>
Foam::label Foam::findLower
(
    const ListType& input,
    const T& val,
    const label start
)
{
    return findLower
    (
        input,
        val,
        start,
        lessOp<T>()
    );
}


template<class ListType>
ListType Foam::reverseList(const ListType& input)
{
    const label len = input.size();
    const label last = (len - 1);

    ListType output(len);
    output.resize(len);     // Consistent sizing (eg, DynamicList)

    for (label i=0; i < len; ++i)
    {
        output[i] = input[last - i];
    }

    return output;
}


template<class ListType>
void Foam::inplaceReverseList(ListType& input)
{
    const label len = input.size();
    const label last = (len - 1);
    const label n2 = len >> 1;

    for (label i=0; i<n2; ++i)
    {
        Foam::Swap(input[i], input[last - i]);
    }
}


template<class ListType>
ListType Foam::rotateList(const ListType& input, const label n)
{
    const label len = input.size();

    ListType output(len);
    output.resize(len);     // Consistent sizing (eg, DynamicList)

    for (label i=0; i<len; ++i)
    {
        label index = (i - n) % len;

        if (index < 0)
        {
            index += len;
        }

        output[i] = input[index];
    }

    return output;
}


template<template<typename> class ListType, class DataType>
void Foam::inplaceRotateList(ListType<DataType>& input, label n)
{
    const label len = input.size();

    n = (len - n) % len;

    if (n < 0)
    {
        n += len;
    }

    SubList<DataType> firstHalf(input, n, 0);
    SubList<DataType> secondHalf(input, len - n, n);

    inplaceReverseList(firstHalf);
    inplaceReverseList(secondHalf);

    inplaceReverseList(input);
}


// * * * * * * * * * * * * * * * * * ListOps * * * * * * * * * * * * * * * * //

template<class T>
void Foam::ListOps::appendEqOp<T>::operator()
(
    List<T>& x,
    const List<T>& y
) const
{
    if (y.size())
    {
        if (x.size())
        {
            x.push_back(y);
        }
        else
        {
            x = y;
        }
    }
}


template<class T>
void Foam::ListOps::uniqueEqOp<T>::operator()
(
    List<T>& x,
    const List<T>& y
) const
{
    if (y.size())
    {
        if (x.size())
        {
            for (const T& val : y)
            {
                x.push_uniq(val);
            }
        }
        else
        {
            x = y;
        }
    }
}


template<class ListType, class UnaryPredicate>
Foam::label Foam::ListOps::count_if
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start
)
{
    label num = 0;

    const label len = input.size();

    if (start >= 0)
    {
        for (label i = start; i < len; ++i)
        {
            if (pred(input[i]))
            {
                ++num;
            }
        }
    }

    return num;
}


template<class ListType, class UnaryPredicate>
Foam::label Foam::ListOps::find_if
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start
)
{
    const label len = input.size();

    if (start >= 0)
    {
        for (label i = start; i < len; ++i)
        {
            if (pred(input[i]))
            {
                return i;
            }
        }
    }

    return -1;
}


template<class ListType, class UnaryPredicate>
bool Foam::ListOps::found_if
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start
)
{
    return (ListOps::find_if(input, pred, start) >= 0);
}


template<class ListType, class UnaryPredicate>
Foam::labelList Foam::ListOps::findIndices
(
    const ListType& input,
    const UnaryPredicate& pred,
    label start
)
{
    const label len = input.size();

    // Pass 1: count occurrences where pred is true. ie, count_if()
    label count = 0;

    if (start >= 0)
    {
        for (label i = start; i < len; ++i)
        {
            if (pred(input[i]))
            {
                if (!count) start = i;  // adjust start for second pass
                ++count;
            }
        }
    }

    labelList indices(count);

    // Pass 2: fill content
    if (count)
    {
        const label total = count;
        count = 0;
        for (label i = start; i < len; ++i)
        {
            if (pred(input[i]))
            {
                indices[count] = i;
                if (++count == total)  // early termination
                {
                    break;
                }
            }
        }
    }

    return indices;
}


template<class T>
void Foam::ListOps::setValue
(
    UList<T>& list,
    const labelUList& locations,
    const T& val
)
{
    const label len = list.size();

    for (const label index : locations)
    {
        // Range-checked
        if (index >= 0 && index < len)
        {
            list[index] = val;
        }
    }
}


template<class T>
void Foam::ListOps::setValue
(
    UList<T>& list,
    const labelHashSet& locations,
    const T& val
)
{
    const label len = list.size();

    for (const label index : locations)
    {
        // Range-checked
        if (index >= 0 && index < len)
        {
            list[index] = val;
        }
    }
}


template<class T>
void Foam::ListOps::setValue
(
    UList<T>& list,
    const UList<bool>& locations,
    const T& val
)
{
    const label len = list.size();
    const label count = locations.size();
    const label end = min(count, len);

    // The efficiency is modest
    for (label index = 0; index < end; ++index)
    {
         if (locations[index])
         {
             list[index] = val;
         }
    }
}


template<class T>
void Foam::ListOps::setValue
(
    UList<T>& list,
    const bitSet& locations,
    const T& val
)
{
    const label len = list.size();

    for
    (
        label pos = locations.find_first();
        pos >= 0 && pos < len;
        pos = locations.find_next(pos)
    )
    {
        list[pos] = val;
    }
}


template<class T, class T2, class UnaryOperation>
Foam::List<T> Foam::ListOps::create
(
    const UList<T2>& input,
    const UnaryOperation& op
)
{
    const label len = input.size();

    List<T> output(len);

    // ie, std::transform(input.begin(), input.end(), output.begin(), op);

    const T2* in = input.begin();
    T* out = output.begin();

    for (label i = 0; i < len; ++i)
    {
        out[i] = op(in[i]);
    }

    return output;
}


template<class T, class InputIterator, class UnaryOperation>
Foam::List<T> Foam::ListOps::create
(
    InputIterator first,
    InputIterator last,
    const UnaryOperation& op
)
{
    const label len = std::distance(first, last);

    List<T> output(len);

    // ie, std::transform(first, last, output.begin(), op);

    T* out = output.begin();

    while (first != last)
    {
        *out = op(*first);
        ++first;
        ++out;
    }

    return output;
}


template<class T>
Foam::List<T> Foam::ListOps::createWithValue
(
    const label len,
    const labelUList& locations,
    const T& val,
    const T& deflt
)
{
    List<T> list(len, deflt);
    ListOps::setValue(list, locations, val);

    return list;
}


template<class T>
Foam::List<T> Foam::ListOps::createWithValue
(
    const label len,
    const labelHashSet& locations,
    const T& val,
    const T& deflt
)
{
    List<T> list(len, deflt);
    ListOps::setValue(list, locations, val);

    return list;
}


template<class T>
Foam::List<T> Foam::ListOps::createWithValue
(
    const label len,
    const UList<bool>& locations,
    const T& val,
    const T& deflt
)
{
    List<T> list(len, deflt);
    ListOps::setValue(list, locations, val);

    return list;
}


template<class T>
Foam::List<T> Foam::ListOps::createWithValue
(
    const label len,
    const bitSet& locations,
    const T& val,
    const T& deflt
)
{
    List<T> list(len, deflt);
    ListOps::setValue(list, locations, val);

    return list;
}


template<class T>
Foam::List<T> Foam::ListOps::createWithValue
(
    const label len,
    const label index,
    const T& val,
    const T& deflt
)
{
    List<T> list(len, deflt);

    // Range-checked
    if (index >= 0 && index < len)
    {
        list[index] = val;
    }

    return list;
}


template<class T>
Foam::List<T> Foam::ListOps::createWithValue
(
    const label len,
    const label index,
    T&& val,
    const T& deflt
)
{
    List<T> list(len, deflt);

    // Range-checked
    if (index >= 0 && index < len)
    {
        list[index] = std::move(val);
    }

    return list;
}


// ************************************************************************* //

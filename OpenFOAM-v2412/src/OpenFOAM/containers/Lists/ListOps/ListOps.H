/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Various functions to operate on Lists.

Namespace
    Foam::ListOps

Description
    Various utility functions to work on Lists.

SourceFiles
    ListOps.C
    ListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ListOps_H
#define Foam_ListOps_H

#include "FlatOutput.H"
#include "labelPair.H"
#include "labelList.H"
#include "IndirectList.H"
#include "HashSet.H"
#include "Map.H"
#include "bitSet.H"
#include "ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class CompactListList;


//- Renumber the values within a list.
//  Negative input values are left untouched.
template<class IntListType>
IntListType renumber(const labelUList& oldToNew, const IntListType& input);

//- Inplace renumber the values within a list.
//  Negative input values are left untouched.
template<class IntListType>
void inplaceRenumber(const labelUList& oldToNew, IntListType& input);

//- Inplace renumber the values of a list.
//  Values that are not mapped by oldToNew are left untouched.
template<class IntListType>
void inplaceRenumber(const Map<label>& oldToNew, IntListType& input);


//- Reorder the elements of a list.
//  Locations with negative oldToNew values are left as is (copy-through).
//  However, if pruning is activated, these negative oldToNew values are
//  instead skipped over and the resulting list shrunk to the max index
//  actually used.
template<class ListType>
ListType reorder
(
    const labelUList& oldToNew,
    const ListType& input,
    const bool prune = false
);

//- Inplace reorder the elements of a list.
//  Locations with negative oldToNew values are left as is (copy-through).
//  However, if pruning is activated, these negative oldToNew values are
//  instead skipped over and the resulting list shrunk to the max index
//  actually used.
template<class ListType>
void inplaceReorder
(
    const labelUList& oldToNew,
    ListType& input,
    const bool prune = false
);


//- Reorder the elements of a packed list.
//  Similar to the general templated form, but with auto-vivify
//  for PackedList.
template<unsigned Width>
PackedList<Width> reorder
(
    const labelUList& oldToNew,
    const PackedList<Width>& input,
    const bool prune = false
);

//- Inplace reorder the elements of a packed list.
//  Similar to the general templated form, but with auto-vivify
//  for PackedList.
template<unsigned Width>
void inplaceReorder
(
    const labelUList& oldToNew,
    PackedList<Width>& input,
    const bool prune = false
);


//- Reorder the elements of a list.
//  Similar to the general templated form, but with auto-vivify for Bitset.
bitSet reorder
(
    const labelUList& oldToNew,
    const bitSet& input,
    const bool prune = false
);

//- Inplace reorder the elements of a list.
//  Similar to the general templated form, but with auto-vivify for bitSet.
void inplaceReorder
(
    const labelUList& oldToNew,
    bitSet& input,
    const bool prune = false
);


//- Rewrite with mapped keys. Ignore elements with negative key.
//  The Container is some type of HashTable or Map with a label for its key.
template<class Container>
void inplaceMapKey(const labelUList& oldToNew, Container& input);

//- Map values. Ignore negative values.
//  \return number of values changed
template<class Container>
label inplaceMapValue(const labelUList& oldToNew, Container& input);

//- Use mapper as a lookup to modify the values of input.
//
//  \return number of values changed
//
//  \code
//  Map<label> mapper({{1, 3}, {2, 6}, {3, 12}, {5, 8}});
//
//  HashTable<label> models with
//  {
//      model0 => 1,
//      model1 => 4,
//      model2 => 5,
//  }
//
//  inplaceMapValue(mapper, models);
//
//  Now contains
//  {
//      model0 => 3,
//      model1 => 4,
//      model2 => 8,
//  }
//  \endcode
//
//  \note the modification occurs in a single pass and will not
//  remap more than once.
template<class Container>
label inplaceMapValue(const Map<label>& mapper, Container& input);


//- Return (sorted) indices corresponding to duplicate list values
template<class T>
labelList duplicateOrder(const UList<T>& input);

//- Generate (sorted) indices corresponding to duplicate list values
template<class T>
void duplicateOrder(const UList<T>& input, labelList& order);

//- Generate (sorted) indices corresponding to duplicate list values
//  sort using specified list compare predicate
template<class T, class ListComparePredicate>
void duplicateOrder
(
    const UList<T>& input,
    labelList& order,
    const ListComparePredicate& comp
);


//- Return (sorted) indices corresponding to unique list values
template<class T>
labelList uniqueOrder(const UList<T>& input);

//- Generate (sorted) indices corresponding to unique list values
template<class T>
void uniqueOrder(const UList<T>& input, labelList& order);

//- Generate (sorted) indices corresponding to unique list values
//  sort using specified list compare predicate
template<class T, class ListComparePredicate>
void uniqueOrder
(
    const UList<T>& input,
    labelList& order,
    const ListComparePredicate& comp
);


//- Return sorted list with removal of duplicates
template<class T>
List<T> uniqueSort(const UList<T>& input);

//- Inplace sorting and removal of duplicates.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class ListType>
void inplaceUniqueSort(ListType& input);

//- Inplace sorting and removal of duplicates.
//  Do not use FixedList for the input list, since it doesn't resize.
template<class ListType, class ListComparePredicate>
void inplaceUniqueSort
(
    ListType& input,
    const ListComparePredicate& comp
);


//- Extract elements of the input list when select is true.
//
//  \param[in] select the bool-list selector, for which the operator[]
//      returns true or false.  A labelHashSet can also be used since
//      it satisfies these requirements
//  \param[in] input the list input values.
//  \param[in] invert set as true to invert the selection logic
//
//  Eg, to extract all selected elements:
//  \code
//    subset<boolList, labelList>(selectedElems, list);
//  \endcode
//
//  \return The subsetted list
//  \see IndirectList::subset
template<class BoolListType, class T>
List<T> subset
(
    const BoolListType& select,
    const UList<T>& input,
    const bool invert=false
);

//- Extract elements of the input list when select is true.
//
//  \param[in] select the selection as a bitSet.
//  \param[in] input the list input values.
//  \param[in] invert set as true to invert the selection logic
//
//  \return The subsetted list
//  \see IndirectList::subset
template<class T>
List<T> subset
(
    const bitSet& select,
    const UList<T>& input,
    const bool invert=false
);

//- Inplace extract elements of the input list when select is true.
//
//  \param[in] select the bool-list selector, for which the operator[]
//      returns true or false.  A labelHashSet can also be used since
//      it satisfies these requirements
//  \param[in,out] input the list input values.
//      Cannot be a FixedList since it doesn't resize.
//  \param[in] invert set as true to invert the selection logic
template<class BoolListType, class ListType>
void inplaceSubset
(
    const BoolListType& select,
    ListType& input,
    const bool invert=false
);

//- Inplace extract elements of the input list when select is true.
//
//  \param[in] select the selection as a bitSet.
//  \param[in,out] input the list input values.
//      Cannot be a FixedList since it doesn't resize.
//  \param[in] invert set as true to invert the selection logic
//
//  \note Includes optimized handling of bitSet when invert=false.
template<class ListType>
void inplaceSubset
(
    const bitSet& select,
    ListType& input,
    const bool invert=false
);

//- Copy a subset of the input list when predicate is true.
//
//  \param[in] input the list input values.
//  \param[in] pred the selection predicate
//  \param[in] invert set as true to invert the selection logic
//
//  \return The subsetted list
//  \see IndirectList::subset_if
template<class T, class UnaryPredicate>
List<T> subsetList
(
    const UList<T>& input,
    const UnaryPredicate& pred,
    const bool invert=false
);


//- Inplace subset of the list when predicate is true.
//
//  \param[in,out] input the list input/output values.
//       Cannot be a FixedList since it doesn't resize.
//  \param[in] pred the selection predicate
//  \param[in] invert set as true to invert the selection logic
template<class ListType, class UnaryPredicate>
void inplaceSubsetList
(
    ListType& input,
    const UnaryPredicate& pred,
    const bool invert=false
);


//- Create an inverse one-to-one mapping.
//  \param len the output size
//  \param map the unique indices to map, in a [0..len) range.
//       Any negative indices are ignored.
//
//  \return the inverse mapping with -1 for unmapped elements.
labelList invert(const label len, const labelUList& map);

//- Create an inverse one-to-one mapping for all 'on' bits of the map.
//  \param len the output size
//  \param map the 'on' bits to use for the mapping indices.
//      The last on bit shall be less than len.
//
//  \return the inverse mapping with -1 for unmapped elements.
labelList invert(const label len, const bitSet& map);

//- Create an inverse one-to-one mapping for all 'on' bits of the map.
//  The output size is dictated by the map size.
//  \param map the unique indices to map.
//
//  \return the inverse mapping with -1 for unmapped elements.
labelList invert(const bitSet& map);

//- Create inverse mapping, which is a lookup table into the given list
//  \param values the values (likely unique).
//
//  Negative values are allowed, duplicates are "first wins"
//
//  \return the inverse lookup
Map<label> invertToMap(const labelUList& values);

//- Invert one-to-many map. Unmapped elements will be size 0.
labelListList invertOneToMany(const label len, const labelUList& map);

//- Invert one-to-many compact map. Unmapped elements will be size 0.
CompactListList<label> invertOneToManyCompact
(
    const label len,
    const labelUList& map
);

//- Invert many-to-many.
//  Input and output types must be inherited from List and also
//  contain ints/labels. Used, for example, for faces to pointFaces.
template<class InputIntListType, class OutputIntListType>
void invertManyToMany
(
    const label len,
    const UList<InputIntListType>& input,
    List<OutputIntListType>& output
);

template<class InputIntListType, class OutputIntListType>
List<OutputIntListType> invertManyToMany
(
    const label len,
    const UList<InputIntListType>& input
)
{
    List<OutputIntListType> output;
    invertManyToMany<InputIntListType,OutputIntListType>(len, input, output);
    return output;
}


//- Deprecated(2017-10) search for first occurrence of the given element.
//  \return The index found or return -1 if not found.
//  \deprecated(2017-10) - use the UList find/found methods
template<class ListType>
FOAM_DEPRECATED_FOR(2017-10, "UList find/found methods")
label findIndex
(
    const ListType& input,
    typename ListType::const_reference val,
    const label start=0
)
{
    return input.find(val, start);
}


//- Linear search to find all occurrences of given element.
template<class ListType>
labelList findIndices
(
    const ListType& input,
    typename ListType::const_reference val,
    label start=0
);


//- Linear search for the index of the min element,
//- similar to std::min_element but for lists and returns the index.
//
//  \tparam ListType  The input list type
//
//  \param input The list to search
//  \param start The start index in the list (default: 0)
//
//  \return The min index or -1 on error.
template<class ListType>
label findMin(const ListType& input, label start=0);

//- Linear search for the index of the max element,
//- similar to std::max_element but for lists and returns the index.
//
//  \tparam ListType  The input list type
//
//  \param input The list to search
//  \param start The start index in the list (default: 0)
//
//  \return The max index or -1 on error.
template<class ListType>
label findMax(const ListType& input, label start=0);


//- Linear search for the index of the min/max element,
//- similar to std::minmax_element but for lists and returns the index.
//
//  \tparam ListType  The input list type
//
//  \param input The list to search
//  \param start The start index in the list (default: 0)
//
//  \return The min/max indices as a Pair (min is first max is second)
//      or (-1,-1) on error.
template<class ListType>
labelPair findMinMax(const ListType& input, label start=0);


//- Binary search to find the index of the last element in a sorted list
//- that is less than value.
//
//  Uses the global <code> < </code> operator and thus
//  <code> (list[i] < val) </code> for the test.
//
//  \tparam ListType  The input list type
//  \tparam T  The value type (should normally be ListType::value_type)
//
//  \param input The sorted list to search
//  \param val   The value for searching/comparing
//  \param start The start index in the list (default: 0)
//
//  \return The index found or -1 if not found.
template<class ListType>
label findSortedIndex
(
    const ListType& input,
    typename ListType::const_reference val,
    const label start=0
);


//- Binary search to find the index of the last element in a sorted list
//- that is less than value.
//
//  Uses <code> lessOp<T>() </code> and thus
//  <code> lessOp<T>(list[i], val) </code> for the test.
//
//  \tparam ListType  The input list type
//  \tparam T  The value type (is often the same as ListType::value_type)
//  \tparam ComparePredicate  The type of the comparison functor that
//      returns true for sorting below.
//
//  \param input The sorted list to search
//  \param val   The value for searching/comparing
//  \param start The start index in the list
//  \param comp  The comparison functor for testing.
//      Uses <code> comp(list[i], val) </code> for the test.
//
//  \return The index found or -1 if not found.
template<class ListType, class T, class ComparePredicate>
label findLower
(
    const ListType& input,
    const T& val,
    const label start,
    const ComparePredicate& comp
);


//- Binary search to find the index of the last element in a sorted list
//- that is less than value.
//
//  Uses <code> lessOp<T>() </code> and thus
//  <code> lessOp<T>(list[i], val) </code> for the test.
//
//  \tparam ListType  The input list type
//  \tparam T  The value type (should normally be ListType::value_type)
//
//  \param input The sorted list to search
//  \param val   The value for searching/comparing
//  \param start The start index in the list (default: 0)
//
//  \return The index found or -1 if not found.
template<class ListType, class T>
label findLower
(
    const ListType& input,
    const T& val,
    const label start=0
);


//- Reverse a list. First element becomes last element etc.
template<class ListType>
ListType reverseList(const ListType& input);


//- Inplace reversal of a list using Swap.
template<class ListType>
void inplaceReverseList(ListType& input);


//- Rotate a list by n places.
//  If n is positive rotate clockwise/right/down.
//  If n is negative rotate anti-clockwise/left/up.
template<class ListType>
ListType rotateList(const ListType& list, const label n);


//- Inplace reversal of a list using the Reversal Block Swapping algorithm.
template<template<typename> class ListType, class DataType>
void inplaceRotateList(ListType<DataType>& list, label n);


/*---------------------------------------------------------------------------*\
                         Namespace ListOps Declaration
\*---------------------------------------------------------------------------*/

namespace ListOps
{

//- List helper to append y elements onto the end of x
template<class T>
struct appendEqOp
{
    void operator()(List<T>& x, const List<T>& y) const;
};

//- List helper to append y unique elements onto the end of x
template<class T>
struct uniqueEqOp
{
    void operator()(List<T>& x, const List<T>& y) const;
};

//- List helper to add y unique elements to x
struct unionEqOp
{
    void operator()(labelList& x, const labelList& y) const;
};


// Public classes

//- A list compare binary predicate for normal sort
template<class ListType>
struct less
{
    const ListType& values;

    less(const ListType& list)
    :
        values(list)
    {}

    bool operator()(const label a, const label b) const
    {
        return (values[a] < values[b]);
    }
};


//- A list compare binary predicate for reverse sort
template<class ListType>
struct greater
{
    const ListType& values;

    greater(const ListType& list)
    :
        values(list)
    {}

    bool operator()(const label a, const label b) const
    {
        return (values[b] < values[a]);
    }
};


//- Fill an identity map with (map[i] == i)
//  Optionally with an alternative start index, so that (map[i] == i+start)
FOAM_DEPRECATED_STRICT(2023-10, "Foam::identity(...)")
inline void identity(labelUList& map, label start = 0)
{
    Foam::identity(map, start);
}


//- Count the number of matching entries.
//  When start is specified, any occurrences before start are ignored.
//  Linear search.
//  Like std::count_if but works with list indexing
template<class ListType, class UnaryPredicate>
label count_if
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start=0
);


//- Find index of the first occurrence that satisfies the predicate.
//  When start is specified, any occurrences before start are ignored.
//  Linear search.
//  Like std::find_if but works with list indexing.
//  \return position in list or -1 if not found.
template<class ListType, class UnaryPredicate>
label find_if
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start=0
);


//- Same as ListOps::find_if
template<class ListType, class UnaryPredicate>
label find
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start=0
)
{
    return ListOps::find_if(input, pred, start);
}


//- True if there is a value in the list that satisfies the predicate.
//  When start is specified, any occurrences before start are ignored.
//  Linear search.
//  \return true if found.
template<class ListType, class UnaryPredicate>
bool found_if
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start=0
);


//- Same as found_if
template<class ListType, class UnaryPredicate>
bool found
(
    const ListType& input,
    const UnaryPredicate& pred,
    const label start=0
)
{
    return ListOps::found_if(input, pred, start);
}


//- Linear search to find all occurences of given element.
template<class ListType, class UnaryPredicate>
labelList findIndices
(
    const ListType& input,
    const UnaryPredicate& pred,
    label start=0
);


//- Set various locations of the list with a specified value.
//
//  \param list the list to modify
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
template<class T>
void setValue
(
    UList<T>& list,
    const labelUList& locations,
    const T& val
);


//- Set various locations of the list with a specified value.
//
//  \param list the list to modify
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
template<class T>
void setValue
(
    UList<T>& list,
    const labelHashSet& locations,
    const T& val
);


//- Set various locations of the list with a specified value.
//
//  \param list the list to modify
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
template<class T>
void setValue
(
    UList<T>& list,
    const UList<bool>& locations,
    const T& val
);


//- Set various locations of the list with a specified value.
//
//  \param list the list to modify
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
template<class T>
void setValue
(
    UList<T>& list,
    const bitSet& locations,
    const T& val
);


//- Create a List from a List of a dissimilar type, using the entire list.
//  For example, convert a list of ints to floats, vectors etc.
//
//  \param input the list input values.
//  \param op the unary conversion operator, which can be used to convert
//      to other types.
//
//  \code
//  auto neg = ListOps::create<label>
//  (
//      ints,
//      std::negate<label>()
//  );
//
//  auto labels = ListOps::create<label>
//  (
//      ints,
//      labelOp<int>()
//  );
//
//  auto vectors = ListOps::create<vector>
//  (
//      ints,
//      [](const int& val){ return vector(1.5*val, 0, 0); }
//  );
//
//  \endcode
template<class T, class T2, class UnaryOperation>
List<T> create
(
    const UList<T2>& input,
    const UnaryOperation& op
);


//- Create a List from an iterator range [first,last) of a dissimilar type.
//  Uses std::distance for the size.
//
//  \param first the begin of the iterator range
//  \param last the end of the iterator range
//  \param op the unary conversion operator, which can be used to convert
//      to other types.
template<class T, class InputIterator, class UnaryOperation>
List<T> create
(
    InputIterator first,
    InputIterator last,
    const UnaryOperation& op
);


//- Create a List filled with default values and various locations with
//- another specified value.
//
//  \param len the length of the list
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
//  \param deflt the initialization default value
template<class T>
List<T> createWithValue
(
    const label len,
    const labelUList& locations,
    const T& val,
    const T& deflt = T()
);


//- Create a List filled with default values and various locations with
//- another specified value.
//
//  \param len the length of the list
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
//  \param deflt the initialization default value
template<class T>
List<T> createWithValue
(
    const label len,
    const labelHashSet& locations,
    const T& val,
    const T& deflt = T()
);


//- Create a List filled with default values and various locations with
//- another specified value.
//
//  \param len the length of the list
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
//  \param deflt the initialization default value
template<class T>
List<T> createWithValue
(
    const label len,
    const UList<bool>& locations,
    const T& val,
    const T& deflt = T()
);


//- Create a List filled with default values and various locations with
//- another specified value.
//
//  \param len the length of the list
//  \param locations where to apply the specified value
//       An out-of-range index is silently ignored.
//  \param val the value to set at the specified locations
//  \param deflt the initialization default value
template<class T>
List<T> createWithValue
(
    const label len,
    const bitSet& locations,
    const T& val,
    const T& deflt = T()
);


//- Create a List filled with default values and one specified value,
//- which is copy assigned at the specified index
//
//  \param len the length of the list
//  \param index where to apply the specified value.
//       An out-of-range index is silently ignored.
//  \param val the value to copy assign at the specified index
//  \param deflt the initialization default value
template<class T>
List<T> createWithValue
(
    const label len,
    const label index,
    const T& val,
    const T& deflt = T()
);


//- Create a List filled with default values and one specified value,
//- which is move assigned at the specified index
//
//  \param len the length of the list
//  \param index where to apply the specified value.
//       An out-of-range index is silently ignored.
//  \param val the value to move assign at the specified index
//  \param deflt the initialization default value
//
//  For example,
//  \code
//  // Gather all unique points on master
//
//  List<pointField> gatheredPoints(Pstream::nProcs());
//  gatheredPoints[Pstream::myProcNo()] = pointField
//  (
//      mesh.points(),
//      uniqueMeshPoints
//  );
//  ...
//
//  // Or else
//  auto gatheredPoints = ListOps::createWithValue<pointField>
//  (
//      Pstream::nProcs(),
//      Pstream::myProcNo(),
//      pointField(mesh.points(), uniqueMeshPoints)
//  );
//  ...
//
//  \endcode
template<class T>
List<T> createWithValue
(
    const label len,
    const label index,
    T&& val,
    const T& deflt = T()
);

} // End namespace ListOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

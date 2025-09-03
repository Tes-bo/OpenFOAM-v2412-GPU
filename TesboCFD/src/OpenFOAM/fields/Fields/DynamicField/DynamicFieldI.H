/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, int SizeMin>
template<class ListType>
inline void Foam::DynamicField<T, SizeMin>::doAssignDynList
(
    const ListType& list
)
{
    const label len = list.size();

    if (capacity_ < len)
    {
        // Needs more space for the copy operation
        List<T>::setAddressableSize(capacity_);  // Use entire space
        List<T>::resize_nocopy(len);
        capacity_ = List<T>::size();
    }

    // Perform copy into addressable portion
    List<T>::setAddressableSize(len);
    List<T>::operator=(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::doCapacity
(
    const bool nocopy,
    const label newCapacity
)
{
    if (newCapacity == capacity_)
    {
        return;
    }

    // Addressable length, possibly truncated by new capacity
    const label currLen = min(List<T>::size(), newCapacity);

    // Corner case - see comments in DynamicList doCapacity
    if (List<T>::size() == newCapacity)
    {
        List<T>::setAddressableSize(currLen+1);
    }

    if (nocopy)
    {
        List<T>::resize_nocopy(newCapacity);
    }
    else
    {
        List<T>::resize(newCapacity);
    }

    capacity_ = List<T>::size();
    List<T>::setAddressableSize(currLen);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::doReserve
(
    const bool nocopy,
    const label len
)
{
    if (capacity_ < len)
    {
        // Preserve addressed size
        const label currLen = List<T>::size();

        // Increase capacity (doubling)
        capacity_ = max(SizeMin, max(len, label(2*capacity_)));

        if (nocopy)
        {
            List<T>::resize_nocopy(capacity_);
        }
        else
        {
            List<T>::resize(capacity_);
        }
        List<T>::setAddressableSize(currLen);
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::doResize
(
    const bool nocopy,
    const label len
)
{
    this->doReserve(nocopy, len);
    List<T>::setAddressableSize(len);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline constexpr Foam::DynamicField<T, SizeMin>::DynamicField() noexcept
:
    Field<T>(),
    capacity_(0)
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField(const label initialCapacity)
:
    Field<T>(),
    capacity_(0)
{
    reserve_nocopy(initialCapacity);
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const label len,
    const T& val
)
:
    Field<T>(len, val),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const label len,
    const Foam::zero
)
:
    Field<T>(len, Foam::zero{}),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const DynamicField<T, SizeMin>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const DynamicField<T, AnySizeMin>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
template<class Addr>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const IndirectListBase<T, Addr>& list
)
:
    Field<T>(list),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    List<T>&& content
) noexcept
:
    Field<T>(std::move(content)),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    DynamicList<T, AnySizeMin>&& content
) noexcept
:
    Field<T>(std::move(static_cast<List<T>&>(content))),
    capacity_(content.capacity())
{
    content.setCapacity_unsafe(0);  // Same as shrink_unsafe() but noexcept
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    DynamicField<T, SizeMin>&& content
) noexcept
:
    Field<T>(std::move(static_cast<List<T>&>(content))),
    capacity_(content.capacity())
{
    content.setCapacity_unsafe(0);  // Same as shrink_unsafe() but noexcept
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    DynamicField<T, AnySizeMin>&& content
) noexcept
:
    Field<T>(std::move(static_cast<List<T>&>(content))),
    capacity_(content.capacity())
{
    content.setCapacity_unsafe(0);  // Same as shrink_unsafe() but noexcept
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& mapF,
    const labelUList& mapAddressing
)
:
    Field<T>(mapF, mapAddressing),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& mapF,
    const labelListList& mapAddressing,
    const scalarListList& weights
)
:
    Field<T>(mapF, mapAddressing, weights),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField
(
    const UList<T>& mapF,
    const FieldMapper& map
)
:
    Field<T>(mapF, map),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>::DynamicField(Istream& is)
:
    Field<T>(is),
    capacity_(Field<T>::size())
{}


template<class T, int SizeMin>
inline Foam::tmp<Foam::DynamicField<T, SizeMin>>
Foam::DynamicField<T, SizeMin>::clone() const
{
    return tmp<DynamicField<T, SizeMin>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline std::streamsize
Foam::DynamicField<T, SizeMin>::capacity_bytes() const noexcept
{
    return std::streamsize(capacity_)*sizeof(T);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::setCapacity
(
    const label len
)
{
    this->doCapacity(false, len);  // nocopy = false
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::setCapacity_nocopy
(
    const label len
)
{
    this->doCapacity(true, len);  // nocopy = true
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::setCapacity_unsafe
(
    const label len
) noexcept
{
    capacity_ = len;
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::reserve
(
    const label len
)
{
    this->doReserve(false, len);  // nocopy = false
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::reserve_nocopy
(
    const label len
)
{
    this->doReserve(true, len);  // nocopy = true
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::resize
(
    const label len
)
{
    this->doResize(false, len);  // nocopy = false
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::resize_nocopy
(
    const label len
)
{
    this->doResize(true, len);  // nocopy = true
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::resize
(
    const label len,
    const T& val
)
{
    const label oldLen = List<T>::size();
    resize(len);

    // Fill newly exposed with constant value
    if (oldLen < List<T>::size())
    {
        std::fill
        (
            this->begin(oldLen), this->end(), val
        );
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::clear() noexcept
{
    List<T>::setAddressableSize(0);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::clearStorage()
{
    List<T>::clear();
    capacity_ = 0;
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::shrink_to_fit()
{
    const label currLen = List<T>::size();

    if (currLen < capacity_)
    {
        // Adjust addressable size to trigger proper resizing
        List<T>::setAddressableSize(currLen+1);
        List<T>::resize(currLen);
        capacity_ = List<T>::size();
    }
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::shrink_unsafe()
{
    if (List<T>::empty())
    {
        // Delete storage if empty
        List<T>::clear();
    }
    capacity_ = List<T>::size();
}


template<class T, int SizeMin>
inline Foam::DynamicField<T, SizeMin>&
Foam::DynamicField<T, SizeMin>::shrink()
{
    this->shrink_to_fit();
    return *this;
}


template<class T, int SizeMin>
inline void
Foam::DynamicField<T, SizeMin>::swap(List<T>& list)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&list)
    )
    {
        return;  // Self-swap is a no-op
    }

    // Remove unused storage
    this->shrink_to_fit();

    // Swap storage and addressable size
    UList<T>::swap(list);

    // Update capacity
    capacity_ = List<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::swap
(
    DynamicField<T, AnySizeMin>& other
)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&other)
    )
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UList<T>::swap(other);

    // Swap capacity
    std::swap(this->capacity_, other.capacity_);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::swap
(
    DynamicList<T, AnySizeMin>& other
)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&other)
    )
    {
        return;  // Self-swap is a no-op
    }

    // Swap storage and addressable size
    UList<T>::swap(other);

    // Swap capacity
    const label oldCap = this->capacity();
    const label newCap = other.capacity();

    this->setCapacity_unsafe(newCap);
    other.setCapacity_unsafe(oldCap);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::transfer(List<T>& list)
{
    Field<T>::transfer(list);
    capacity_ = Field<T>::size();
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::transfer
(
    DynamicList<T, AnySizeMin>& list
)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&list)
    )
    {
        return;  // Self-assignment is a no-op
    }

    // Take over storage as-is (without shrink)
    capacity_ = list.capacity();
    Field<T>::transfer(static_cast<List<T>&>(list));
    list.clearStorage();  // capacity=0 etc.
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::transfer
(
    DynamicField<T, AnySizeMin>& list
)
{
    if
    (
        static_cast<const List<T>*>(this)
     == static_cast<const List<T>*>(&list)
    )
    {
        return;  // Self-assignment is a no-op
    }

    // Take over storage as-is (without shrink)
    capacity_ = list.capacity();
    Field<T>::transfer(static_cast<List<T>&>(list));
    list.clearStorage();  // capacity=0 etc.
}


template<class T, int SizeMin>
template<class... Args>
inline T& Foam::DynamicField<T, SizeMin>::emplace_back(Args&&... args)
{
    // This could/should be better with inplace construction
    // (as per std::vector), but currently lacking the methods for that
    // so resize and move assign

    const label idx = List<T>::size();
    resize(idx + 1);

    // move assign element
    this->operator[](idx) = T(std::forward<Args>(args)...);
    return this->operator[](idx);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::push_back
(
    const T& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    this->operator[](idx) = val;  // copy element
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::push_back
(
    T&& val
)
{
    const label idx = List<T>::size();
    resize(idx + 1);

    this->operator[](idx) = std::move(val);  // move assign element
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::push_back
(
    const UList<T>& list
)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    std::copy(list.begin(), list.end(), this->begin(idx));
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::push_back
(
    List<T>&& list
)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted push_back to self"
            << abort(FatalError);
    }

    const label idx = List<T>::size();
    resize(idx + list.size());

    std::move(list.begin(), list.end(), this->begin(idx));

    list.clear();
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::pop_back(label n)
{
    if (n >= this->size())
    {
        this->clear();
    }
    else if (n > 0)
    {
        resize(this->size() - n);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, int SizeMin>
inline T& Foam::DynamicField<T, SizeMin>::operator()
(
    const label i
)
{
    if (i >= List<T>::size())
    {
        resize(i + 1);
    }

    return this->operator[](i);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const T& val
)
{
    UList<T>::operator=(val);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const Foam::zero
)
{
    UList<T>::operator=(Foam::zero{});
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const UList<T>& list
)
{
    doAssignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const DynamicField<T, SizeMin>& list
)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    doAssignDynList(list);
}


template<class T, int SizeMin>
template<class Addr>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    const IndirectListBase<T, Addr>& list
)
{
    // NOTE: Self-assignment needs special handling
    /// if
    /// (
    ///     static_cast<const UList<T>*>(this)
    ///  == static_cast<const UList<T>*>(&list.values())
    /// )

    doAssignDynList(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    List<T>&& list
)
{
    transfer(list);
}


template<class T, int SizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    DynamicField<T, SizeMin>&& list
)
{
    transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    DynamicField<T, AnySizeMin>&& list
)
{
    transfer(list);
}


template<class T, int SizeMin>
template<int AnySizeMin>
inline void Foam::DynamicField<T, SizeMin>::operator=
(
    DynamicList<T, AnySizeMin>&& list
)
{
    transfer(list);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, int SizeMin>
inline Foam::Istream& Foam::DynamicField<T, SizeMin>::readList
(
    Istream& is
)
{
    // Use DynamicList::readList for reading DynamicField.
    // The logic should be the same and this avoids duplicate code

    DynamicList<T, SizeMin> list;
    this->swap(list);

    list.readList(is);
    this->swap(list);

    return is;
}


template<class T, int SizeMin>
inline Foam::Istream& Foam::operator>>
(
    Istream& is,
    DynamicField<T, SizeMin>& rhs
)
{
    return rhs.readList(is);
}


template<class T, int SizeMin>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DynamicField<T, SizeMin>& rhs
)
{
    os << static_cast<const Field<T>&>(rhs);
    return os;
}


// ************************************************************************* //

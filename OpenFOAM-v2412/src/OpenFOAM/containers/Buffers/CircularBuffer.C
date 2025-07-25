/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

template<class T>
void Foam::CircularBuffer<T>::doReserve
(
    const bool nocopy,
    const label len
)
{
    if (storage_.size() < len)
    {
        // Increase capacity (doubling)
        const label newCapacity =
            max(min_size(), max(len+1, label(2*storage_.size())));

        if (nocopy || empty())
        {
            // Simple - no content to preserve

            clear();  // Reset begin/end
            storage_.resize_nocopy(newCapacity);
        }
        else
        {
            // Preserve content
            const labelRange range1 = range_one();
            const labelRange range2 = range_two();

            List<T> old(newCapacity);
            storage_.swap(old);
            begin_ = 0;
            end_ = 0;

            for (const label i : range1)
            {
                storage_[end_++] = std::move(old[i]);
            }
            for (const label i : range2)
            {
                storage_[end_++] = std::move(old[i]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::SubList<T> Foam::CircularBuffer<T>::array_one()
{
    const label len = size_one();
    return (len ? storage_.slice(begin_, len) : SubList<T>());
}


template<class T>
Foam::SubList<T> Foam::CircularBuffer<T>::array_two()
{
    const label len = size_two();
    return (len ? storage_.slice(0, len) : SubList<T>());
}


template<class T>
const Foam::SubList<T> Foam::CircularBuffer<T>::array_one() const
{
    const label len = size_one();
    return (len ? storage_.slice(begin_, len) : SubList<T>());
}


template<class T>
const Foam::SubList<T> Foam::CircularBuffer<T>::array_two() const
{
    const label len = size_two();
    return (len ? storage_.slice(0, len) : SubList<T>());
}


template<class T>
Foam::label Foam::CircularBuffer<T>::find(const T& val, label pos) const
{
    if (pos < 0) return -1;  // no-op

    label i = -1;

    const auto list1 = this->array_one();

    if (pos < list1.size())
    {
        // Can start search in first array
        i = list1.find(val, pos);

        // Position for continued search in second array
        pos = 0;
    }
    else
    {
        // Position for continued search in second array
        pos -= list1.size();
    }

    const auto list2 = this->array_two();

    if (i < 0 && list2.size())
    {
        // Not yet found, continue search in second array
        i = list2.find(val, pos);

        if (i >= 0)
        {
            // As flat index into the entire buffer
            i += list1.size();
        }
    }

    return i;
}


template<class T>
void Foam::CircularBuffer<T>::reverse()
{
    const label n = this->size();
    const label nBy2 = n/2;

    for (label i = 0; i < nBy2; ++i)
    {
        Foam::Swap(operator[](i), operator[](n-1-i));
    }
}


template<class T>
Foam::List<T> Foam::CircularBuffer<T>::list() const
{
    const auto list1 = array_one();
    const auto list2 = array_two();

    List<T> result(list1.size() + list2.size());

    if (list1.size())
    {
        result.slice(0, list1.size()) = list1;
    }
    if (list2.size())
    {
        result.slice(list1.size(), list1.size() + list2.size()) = list2;
    }

    return result;
}


// ************************************************************************* //

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

Class
    Foam::LPtrList

Description
    Template class for non-intrusive linked PtrLists.

SourceFiles
    LPtrList.C
    LPtrListIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LPtrList_H
#define Foam_LPtrList_H

#include "LList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class LListBase, class T> class LPtrList;

template<class LListBase, class T>
Istream& operator>>
(
    Istream& is,
    LPtrList<LListBase, T>& list
);

template<class LListBase, class T>
Ostream& operator<<
(
    Ostream& os,
    const LPtrList<LListBase, T>& list
);


/*---------------------------------------------------------------------------*\
                          Class LPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class LPtrList
:
    public LList<LListBase, T*>
{
    // Private Member Functions

        //- Read from Istream using given Istream constructor class
        template<class INew>
        void readIstream(Istream& is, const INew& inew);


public:

    // STL type definitions

        //- Pointer for LPtrList::value_type objects.
        typedef T* pointer;

        //- Const pointer for LPtrList::value_type objects.
        typedef const T* const_pointer;

        //- Reference for LPtrList::value_type objects.
        typedef T& reference;

        //- Const reference for LPtrList::value_type objects.
        typedef const T& const_reference;


    // Forward Declaration (iterators)

        class iterator;
        class const_iterator;

        using base_iterator = typename LListBase::iterator;
        using const_base_iterator = typename LListBase::const_iterator;

        //- The parent list storage
        typedef LList<LListBase, T*> parent_type;


    // Constructors

        //- Default construct
        LPtrList() = default;

        //- Construct and add initial item pointer
        explicit LPtrList(T* item)
        {
            this->push_front(item);
        }

        //- Copy construct by using 'clone()' for each element
        LPtrList(const LPtrList& lst);

        //- Move construct
        LPtrList(LPtrList&& lst);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        LPtrList(Istream& is, const INew& inew);

        //- Construct from Istream using default Istream constructor class
        explicit LPtrList(Istream& is);


    //- Destructor. Calls clear()
    ~LPtrList();


    // Member Functions

        //- The first entry in the list
        T& front()
        {
            return *(parent_type::front());
        }

        //- The first entry in the list (const access)
        const T& front() const
        {
            return *(parent_type::front());
        }

        //- The last entry in the list
        T& back()
        {
            return *(parent_type::back());
        }

        //- The last entry in the list (const access)
        const T& back() const
        {
            return *(parent_type::back());
        }

        //- Emplace construct new element at front of list. Return reference.
        template<class... Args>
        T& emplace_front(Args&&... args)
        {
            T* ptr = new T(std::forward<Args>(args)...);
            this->push_front(ptr);
            return *ptr;  // OR: return *(parent_type::front());
        }

        //- Emplace construct new element at back of list. Return reference.
        template<class... Args>
        T& emplace_back(Args&&... args)
        {
            T* ptr = new T(std::forward<Args>(args)...);
            this->push_back(ptr);
            return *ptr;  // OR: return *(parent_type::back());
        }

        //- Remove first element(s) from the list (deletes pointers)
        void pop_front(label n = 1);

        //- Clear the contents of the list
        void clear();

        //- Transfer the contents of the argument into this List
        //- and annul the argument list.
        void transfer(LPtrList<LListBase, T>& lst);


    // Member Operators

        //- Copy assign by using 'clone()' for each element
        void operator=(const LPtrList<LListBase, T>& lst);

        //- Move assign
        void operator=(LPtrList<LListBase, T>&& lst);


    // STL iterator

        //- An STL-conforming iterator
        class iterator
        :
            public parent_type::iterator
        {
        public:

            iterator(base_iterator iter)
            :
                parent_type::iterator(iter)
            {}

            //- Return the address of the object being referenced
            pointer get() const
            {
                return parent_type::iterator::operator*();
            }

            reference operator*() const
            {
                return *(this->get());
            }

            pointer operator->() const
            {
                return this->get();
            }

            reference operator()() const
            {
                return operator*();
            }
        };


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public parent_type::const_iterator
        {
        public:

            const_iterator(const_base_iterator iter)
            :
                parent_type::const_iterator(iter)
            {}

            const_iterator(base_iterator iter)
            :
                parent_type::const_iterator(iter)
            {}

            //- Return the address of the object being referenced
            const_pointer get() const
            {
                return parent_type::const_iterator::operator*();
            }

            const_reference operator*() const
            {
                return *(this->get());
            }

            const_pointer operator->() const
            {
                return this->get();
            }

            const_reference operator()() const
            {
                return operator*();
            }
        };


    // STL reverse_iterator

        //- A reverse_iterator, for base classes that support
        //- reverse iteration
        class reverse_iterator
        :
            public parent_type::reverse_iterator
        {
        public:

            reverse_iterator(base_iterator iter)
            :
                parent_type::reverse_iterator(iter)
            {}

            //- Return the address of the object being referenced
            pointer get() const
            {
                return parent_type::reverse_iterator::operator*();
            }

            reference operator*() const
            {
                return *(this->get());
            }

            pointer operator->() const
            {
                return this->get();
            }

            reference operator()() const
            {
                return operator*();
            }
        };


    // STL const_reverse_iterator

        //- A const_reverse_iterator, for base classes that support
        //- reverse iteration
        class const_reverse_iterator
        :
            public parent_type::const_reverse_iterator
        {
        public:

            const_reverse_iterator(const_base_iterator iter)
            :
                parent_type::const_reverse_iterator(iter)
            {}

            //- Return the address of the object being referenced
            const_pointer get() const
            {
                return parent_type::const_reverse_iterator::operator*();
            }

            const_reference operator*() const
            {
                return *(this->get());
            }

            const_pointer operator->() const
            {
                return this->get();
            }

            const_reference operator()() const
            {
                return operator*();
            }
        };


        //- Iterator to first item in list with non-const access
        inline iterator begin()
        {
            return LListBase::template iterator_first<base_iterator>();
        }

        //- Iterator to first item in list with const access
        inline const_iterator cbegin() const
        {
            return LListBase::template iterator_first<const_base_iterator>();
        }

        //- Iterator to last item in list with non-const access
        inline reverse_iterator rbegin()
        {
            return LListBase::template iterator_last<base_iterator>();
        }

        //- Iterator to last item in list with const access
        inline const_reverse_iterator crbegin() const
        {
            return LListBase::template iterator_last<const_base_iterator>();
        }

        //- Iterator to first item in list with const access
        inline const_iterator begin() const
        {
            return LListBase::cbegin();
        }

        //- Iterator to last item in list with const access
        inline const_reverse_iterator rbegin() const
        {
            return crbegin();
        }


        //- End of list for forward iterators
        inline const iterator& end()
        {
            return LListBase::template iterator_end<iterator>();
        }

        //- End of list for forward iterators
        inline const const_iterator& cend() const
        {
            return LListBase::template iterator_end<const_iterator>();
        }

        //- End of list for reverse iterators
        inline const reverse_iterator& rend()
        {
            return LListBase::template iterator_rend<reverse_iterator>();
        }

        //- End of list for reverse iterators
        inline const const_reverse_iterator& crend() const
        {
            return LListBase::template iterator_rend<const_reverse_iterator>();
        }

        //- End of list for forward iterators
        inline const const_iterator& end() const
        {
            return cend();
        }

        //- End of list for reverse iterators
        inline const const_reverse_iterator& rend() const
        {
            return crend();
        }


    // IOstream operators

        friend Istream& operator>> <LListBase, T>
        (
            Istream& is,
            LPtrList<LListBase, T>& list
        );

        friend Ostream& operator<< <LListBase, T>
        (
            Ostream& os,
            const LPtrList<LListBase, T>& list
        );


    // Housekeeping

        //- The first entry in the list
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        T& first() { return front(); }

        //- The first entry in the list (const access)
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        const T& first() const { return front(); }

        //- The last entry in the list
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        T& last() { return back(); }

        //- The last entry in the list (const access)
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        const T& last() const { return back(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LPtrList.C"
    #include "LPtrListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::UILList

Description
    Template class for intrusive linked lists.

SourceFiles
    UILList.C
    UILListIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_UILList_H
#define Foam_UILList_H

#include "label.H"
#include "uLabel.H"
#include "stdFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class LListBase, class T> class UILList;

template<class LListBase, class T>
Ostream& operator<<
(
    Ostream& os,
    const UILList<LListBase, T>& lst
);


/*---------------------------------------------------------------------------*\
                           Class UILList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class UILList
:
    public LListBase
{
public:

    // STL type definitions

        //- Type of values stored
        typedef T value_type;

        //- Pointer for value_type
        typedef T* pointer;

        //- Const pointer for value_type
        typedef const T* const_pointer;

        //- Reference for value_type
        typedef T& reference;

        //- Const reference for value_type
        typedef const T& const_reference;

        //- The type that can represent the container size
        typedef label size_type;

        //- The difference between iterator objects
        typedef label difference_type;


    // Forward Declarations (iterators)

        class iterator;
        class const_iterator;

        using base_iterator = typename LListBase::iterator;
        using const_base_iterator = typename LListBase::const_iterator;


    // Constructors

        //- Default construct
        UILList() = default;

        //- Construct and add initial item pointer
        explicit UILList(T* item)
        {
            this->push_front(item);
        }

        //- Construct as copy
        UILList(const UILList<LListBase, T>& list);


    // Member Functions

        //- The first entry in the list
        T* front()
        {
            return static_cast<T*>(LListBase::front());
        }

        //- The first entry in the list (const access)
        const T* front() const
        {
            return static_cast<const T*>(LListBase::front());
        }

        //- The last entry in the list
        T* back()
        {
            return static_cast<T*>(LListBase::back());
        }

        //- The last entry in the list (const access)
        const T* back() const
        {
            return static_cast<const T*>(LListBase::back());
        }


        //- Remove and return head
        T* removeHead()
        {
            return static_cast<T*>(LListBase::removeHead());
        }

        //- Remove and return element
        T* remove(T* item)
        {
            return static_cast<T*>(LListBase::remove(item));
        }

        //- Remove and return item specified by iterator
        T* remove(iterator& iter)
        {
            return static_cast<T*>(LListBase::remove(iter));
        }


    // Member Operators

        //- Copy assignment
        void operator=(const UILList<LListBase, T>& lst);

        //- Equality. True both lists are element-wise equal
        //  (using value_type::operator==).  Takes linear time.
        bool operator==(const UILList<LListBase, T>& lst) const;

        //- The opposite of the equality operation. Takes linear time.
        bool operator!=(const UILList<LListBase, T>& lst) const;


    // IOstream Operators

        //- Write UILList with line-breaks when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;

        //- Write UILList to Ostream with line breaks,
        //- as per writeList() with shortLen=-1
        friend Ostream& operator<< <LListBase, T>
        (
            Ostream& os,
            const UILList<LListBase, T>& lst
        );


    // STL iterator

        //- A non-const iterator
        class iterator
        :
            public base_iterator
        {
        public:

            iterator(base_iterator iter)
            :
                base_iterator(iter)
            {}

            //- Return the address of the object being referenced
            pointer get() const
            {
                return static_cast<T*>(base_iterator::get_node());
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

            iterator& operator++()
            {
                this->next();
                return *this;
            }
        };


    // STL const_iterator

        //- A const_iterator
        class const_iterator
        :
            public const_base_iterator
        {
        public:

            //- Construct from base const_iterator
            const_iterator(const_base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            //- Construct from base iterator
            const_iterator(base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            //- Return the address of the object being referenced
            const_pointer get() const
            {
                return static_cast<const T*>(const_base_iterator::get_node());
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

            const_iterator& operator++()
            {
                this->next();
                return *this;
            }
        };


    // STL reverse_iterator

        //- A reverse_iterator, for LListBase classes that support
        //- reverse iteration
        class reverse_iterator
        :
            public base_iterator
        {
        public:

            reverse_iterator(base_iterator iter)
            :
                base_iterator(iter)
            {}

            //- Return the address of the object being referenced
            pointer get() const
            {
                return static_cast<T*>(base_iterator::get_node());
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

            reverse_iterator& operator++()
            {
                this->prev();   // Only if base iterator is bidirectional
                return *this;
            }
        };


    // STL const_reverse_iterator

        //- A const_reverse_iterator, for LListBase classes that support
        //- reverse iteration
        class const_reverse_iterator
        :
            public const_base_iterator
        {
        public:

            const_reverse_iterator(const_base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            //- Return the address of the object being referenced
            const_pointer get() const
            {
                return static_cast<const T*>(const_base_iterator::get_node());
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

            const_reverse_iterator& operator++()
            {
                this->prev();   // Only if base iterator is bidirectional
                return *this;
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


    // Housekeeping

        //- The first entry in the list
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        T* first() { return front(); }

        //- The first entry in the list (const access)
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        const T* first() const { return front(); }

        //- The last entry in the list
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        T* last() { return back(); }

        //- The last entry in the list (const access)
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        const T* last() const { return back(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UILList.C"
    #include "UILListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

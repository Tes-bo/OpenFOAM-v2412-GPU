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
    Foam::LList

Description
    Template class for non-intrusive linked lists.

SourceFiles
    LList.C
    LListIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_LList_H
#define Foam_LList_H

#include "label.H"
#include "stdFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class LListBase, class T> class LList;

template<class LListBase, class T>
Istream& operator>>
(
    Istream& is,
    LList<LListBase, T>& lst
);

template<class LListBase, class T>
Ostream& operator<<
(
    Ostream& os,
    const LList<LListBase, T>&  lst
);


/*---------------------------------------------------------------------------*\
                            Class LList Declaration
\*---------------------------------------------------------------------------*/

template<class LListBase, class T>
class LList
:
    public LListBase
{
public:

    // STL type definitions

        //- Type of values stored.
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

        //- The difference between iterators
        typedef label difference_type;


    // Forward Declarations (iterators)

        class iterator;
        class const_iterator;

        using base_iterator = typename LListBase::iterator;
        using const_base_iterator = typename LListBase::const_iterator;


    //- The storage of T with linked nodes
    struct link
    :
        public LListBase::link
    {
        //- Stored object
        T val_;

        //- Copy construct from given object
        link(const T& elem)
        :
            val_(elem)
        {}

        //- Move construct from given object
        link(T&& elem)
        :
            val_(std::move(elem))
        {}


        //- Delete linked item and return the element value
        static T remove(typename LListBase::link* node)
        {
            link* p = static_cast<link*>(node);
            T val(std::move(p->val_));
            delete p;
            return val;
        }

        //- Dereference LListBase::link to obtain address of stored object
        static constexpr T* ptr(typename LListBase::link* node)
        {
            return &(static_cast<link*>(node)->val_);
        }

        //- Dereference LListBase::link to obtain address of stored object
        static constexpr const T* ptr(const typename LListBase::link* node)
        {
            return &(static_cast<const link*>(node)->val_);
        }

        //- Dereference LListBase::link to obtain the stored object
        static constexpr T& ref(typename LListBase::link* node)
        {
            return static_cast<link*>(node)->val_;
        }

        //- Dereference LListBase::link to obtain the stored object
        static constexpr const T& ref(const typename LListBase::link* node)
        {
            return static_cast<const link*>(node)->val_;
        }
    };


    // Constructors

        //- Default construct
        LList() = default;

        //- Construct and copy add initial item
        explicit LList(const T& elem)
        {
            this->push_front(elem);
        }

        //- Construct and move add initial item
        explicit LList(T&& elem)
        {
            this->push_front(std::move(elem));
        }

        //- Construct from Istream
        explicit LList(Istream& is);

        //- Copy construct
        LList(const LList<LListBase, T>& lst);

        //- Move construct
        LList(LList<LListBase, T>&& lst);

        //- Copy construct from an initializer list
        LList(std::initializer_list<T> lst);


    //- Destructor. Calls clear()
    ~LList();


    // Member Functions

        //- The first entry in the list
        reference front()
        {
            return link::ref(LListBase::front());
        }

        //- The first entry in the list (const access)
        const_reference front() const
        {
            return link::ref(LListBase::front());
        }

        //- The last entry in the list
        reference back()
        {
            return link::ref(LListBase::back());
        }

        //- The last entry in the list (const access)
        const_reference back() const
        {
            return link::ref(LListBase::back());
        }


        //- Add copy at front of list
        void push_front(const T& elem)
        {
            LListBase::push_front(new link(elem));
        }

        //- Move construct at front of list
        void push_front(T&& elem)
        {
            LListBase::push_front(new link(std::move(elem)));
        }

        //- Add copy at back of list
        void push_back(const T& elem)
        {
            LListBase::push_back(new link(elem));
        }

        //- Move construct at back of list
        void push_back(T&& elem)
        {
            LListBase::push_back(new link(std::move(elem)));
        }

        //- Delete contents of list
        void clear();

        //- Remove first element(s) from the list (deletes pointers)
        void pop_front(label n = 1);

        //- Remove and return first entry
        T removeHead()
        {
            return link::remove(LListBase::removeHead());
        }

        //- Remove and return element
        T remove(link* item)
        {
            return link::remove(LListBase::remove(item));
        }

        //- Remove and return element specified by iterator
        T remove(iterator& iter)
        {
            return link::remove(LListBase::remove(iter));
        }

        //- Transfer the contents of the argument into this List
        //- and annul the argument list.
        void transfer(LList<LListBase, T>& lst);


    // Member Operators

        //- Copy assignment
        void operator=(const LList<LListBase, T>& lst);

        //- Move assignment
        void operator=(LList<LListBase, T>&& lst);

        //- Copy assignment from an initializer list
        void operator=(std::initializer_list<T> lst);


    // IOstream Operators

        //- Read list from Istream
        Istream& readList(Istream& is);

        //- Write LList with line-breaks when length exceeds shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeList(Ostream& os, const label shortLen=0) const;

        //- Read list from Istream
        friend Istream& operator>> <LListBase, T>
        (
            Istream&,
            LList<LListBase, T>& lst
        );

        //- Write LList to Ostream with line breaks,
        //- as per writeList with shortLen=-1
        friend Ostream& operator<< <LListBase, T>
        (
            Ostream& os,
            const LList<LListBase, T>& lst
        );


    // STL iterator

        //- An STL-conforming iterator
        class iterator
        :
            public base_iterator
        {
        public:

            //- Construct from base iterator
            iterator(base_iterator iter)
            :
                base_iterator(iter)
            {}

            reference operator*() const
            {
                return link::ref(this->get_node());
            }

            pointer operator->() const
            {
                return link::ptr(this->get_node());
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

            iterator& operator--()
            {
                this->prev();   // May not be implemented
                return *this;
            }
        };


    // STL const_iterator

        //- An STL-conforming const_iterator
        class const_iterator
        :
            public const_base_iterator
        {
        public:

            //- Construct from base iterator
            const_iterator(const_base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            //- Construct from base iterator
            const_iterator(base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            const_reference operator*() const
            {
                return link::ref(this->get_node());
            }

            const_pointer operator->() const
            {
                return link::ptr(this->get_node());
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

            const_iterator& operator--()
            {
                this->prev();   // May not be implemented
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

            //- Construct from base iterator
            reverse_iterator(base_iterator iter)
            :
                base_iterator(iter)
            {}

            reference operator*() const
            {
                return link::ref(this->get_node());
            }

            pointer operator->() const
            {
                return link::ptr(this->get_node());
            }

            reverse_iterator& operator++()
            {
                this->prev();   // Only if base iterator is bidirectional
                return *this;
            }

            reverse_iterator& operator--()
            {
                this->next();
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

            //- Construct from base iterator
            const_reverse_iterator(const_base_iterator iter)
            :
                const_base_iterator(iter)
            {}

            const_reference operator*() const
            {
                return link::ref(this->get_node());
            }

            const_pointer operator->() const
            {
                return link::ptr(this->get_node());
            }

            const_reverse_iterator& operator++()
            {
                this->prev();   // Only if base iterator is bidirectional
                return *this;
            }

            const_reverse_iterator& operator--()
            {
                this->next();
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
        reference first() { return front(); }

        //- The first entry in the list (const access)
        //FOAM_DEPRECATED_FOR(2022-10, "front()")
        const_reference first() const { return front(); }

        //- The last entry in the list
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        reference last() { return back(); }

        //- The last entry in the list (const access)
        //FOAM_DEPRECATED_FOR(2022-10, "back()")
        const_reference last() const { return back(); }

        //- Add copy at front of list
        //FOAM_DEPRECATED_FOR(2022-10, "push_front()")
        void prepend(const T& elem) { push_front(elem); }

        //- Move construct at front of list
        //FOAM_DEPRECATED_FOR(2022-10, "push_front()")
        void prepend(T&& elem) { push_front(std::move(elem)); }

        //- Add copy at back of list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(const T& elem) { push_back(elem); }

        //- Move construct at back of list
        //FOAM_DEPRECATED_FOR(2022-10, "push_back()")
        void append(T&& elem) { push_back(std::move(elem)); }

        //- Add copy at front of list. Same as push_front()
        //FOAM_DEPRECATED_FOR(2022-01, "push_front()")
        void insert(const T& elem) { push_front(elem); }

        //- Move construct at front of list. Same as push_front()
        //FOAM_DEPRECATED_FOR(2022-01, "push_front()")
        void insert(T&& elem) { push_front(std::move(elem)); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LList.C"
    #include "LListIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

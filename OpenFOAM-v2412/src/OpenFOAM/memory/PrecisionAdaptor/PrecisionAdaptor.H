/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::PrecisionAdaptor

Description
    Conversion adaptor for Field/List that either wrap the input as a
    reference, or creates a temporary pointer and copies the values
    on construction/destruction.

    This provides, for example, automatic conversion between types
    for linear solvers able to run mixed precision.

\*---------------------------------------------------------------------------*/

#ifndef Foam_PrecisionAdaptor_H
#define Foam_PrecisionAdaptor_H

#include <algorithm>    // For std::copy
#include <type_traits>  // For std::is_same
#include "refPtr.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ConstPrecisionAdaptor Declaration
\*---------------------------------------------------------------------------*/

//- A const Field/List wrapper with possible data conversion
template<class Type, class InputType, template<class> class Container = Field>
class ConstPrecisionAdaptor
:
    public refPtr<Container<Type>>
{
    // Private Member Functions

        //- Set adaptor for different input, copying as required
        void setInput(const Container<InputType>& src)
        {
            if (std::is_same<Type, InputType>::value)
            {
                // Use reference directly
                this->cref(reinterpret_cast<const Container<Type>&>(src));
            }
            else
            {
                // Need intermediate buffer
                this->reset(new Container<Type>(src.size()));
                std::copy(src.cbegin(), src.cend(), this->ref().begin());
            }
        }

        //- Set from tmp, steal pointer if possible
        void tmpInput(tmp<Container<InputType>>& tsrc)
        {
            if (std::is_same<Type, InputType>::value && tsrc.is_pointer())
            {
                // Acquire control of the managed pointer
                this->reset(reinterpret_cast<Container<Type>*>(tsrc.ptr()));
            }
            else
            {
                this->setInput(tsrc.cref());
            }
            tsrc.clear();
        }

public:

    //- The adapted field type. Same as element_type
    typedef Container<Type> FieldType;


    // Constructors

        //- Default construct, setting content later
        ConstPrecisionAdaptor() = default;

        //- Construct from Container of InputType, copying if required
        explicit ConstPrecisionAdaptor(const Container<InputType>& input)
        {
            this->setInput(input);
        }

        //- Construct from tmp Container of InputType, copy/move as required
        explicit ConstPrecisionAdaptor(tmp<Container<InputType>>&& input)
        {
            this->tmpInput(input);
        }

        //- Construct from tmp Container of InputType, copy/move as required
        explicit ConstPrecisionAdaptor(const tmp<Container<InputType>>& input)
        {
            this->tmpInput(const_cast<tmp<Container<InputType>>&>(input));
        }


    // Member Functions

        //- Is precision adaption being used (non-passive adaptor)?
        bool active() const noexcept
        {
            // Same as refPtr::movable()
            return (this->is_pointer() && this->good());
        }

        //- Commit adapted content changes (no-op for const adaptor)
        void commit()
        {}

        //- Set adaptor for different input, copying input if required
        void set(const Container<InputType>& input)
        {
            this->setInput(input);
        }

        //- Set adaptor for tmp Container of InputType, copy/move as required
        void set(tmp<Container<InputType>>&& input)
        {
            this->tmpInput(input);
        }

        //- Set adaptor for tmp Container of InputType, copy/move as required
        void set(const tmp<Container<InputType>>& input)
        {
            this->tmpInput(const_cast<tmp<Container<InputType>>&>(input));
        }


    // Static Member Functions

        //- Select a reference to the input (if types are identical),
        //- or copy into other and return a reference to that
        static const Container<Type>& select
        (
            const Container<InputType>& input,
            Container<Type>& other
        )
        {
            if (std::is_same<Type, InputType>::value)
            {
                return reinterpret_cast<const Container<Type>&>(input);
            }
            else
            {
                other.resize(input.size());
                std::copy(input.cbegin(), input.cend(), other.begin());
                return other;
            }
        }
};


/*---------------------------------------------------------------------------*\
                      Class PrecisionAdaptor Declaration
\*---------------------------------------------------------------------------*/

//- A non-const Field/List wrapper with possible data conversion
template<class Type, class InputType, template<class> class Container = Field>
class PrecisionAdaptor
:
    public refPtr<Container<Type>>
{
    // Private Data

        //- Reference to underlying external input data
        refPtr<Container<InputType>> orig_;


    // Private Member Functions

        //- Set adaptor for different input, copying as required
        void setInput(Container<InputType>& src, const bool doCopy)
        {
            orig_.ref(src);
            if (std::is_same<Type, InputType>::value)
            {
                // Use reference directly
                this->ref(reinterpret_cast<Container<Type>&>(src));
            }
            else
            {
                // Need intermediate buffer
                this->reset(new Container<Type>(src.size()));
                if (doCopy)
                {
                    std::copy(src.cbegin(), src.cend(), this->ref().begin());
                }
            }
        }

public:

    //- The adapted field type. Same as element_type
    typedef Container<Type> FieldType;


    // Constructors

        //- Default construct, setting content later
        PrecisionAdaptor() = default;

        //- Construct from Container<InputType>,
        //- copying input if required (and requested)
        explicit PrecisionAdaptor
        (
            Container<InputType>& input,
            const bool doCopy = true
        )
        {
            this->setInput(input, doCopy);
        }


    //- Destructor, copies back content changes (as required)
    ~PrecisionAdaptor()
    {
        this->commit();  // Commit changes
        this->clear();
    }


    // Member Functions

        //- Is precision adaption being used (non-passive adaptor)?
        bool active() const noexcept
        {
            // Same as refPtr::movable()
            return (this->is_pointer() && this->good());
        }

        //- Commit adapted content changes back to original input (as required)
        void commit()
        {
            if (this->active() && orig_.good())
            {
                const auto& stored = this->cref();
                auto& input = orig_.ref();
                input.resize(stored.size());  // Extra safety
                std::copy(stored.cbegin(), stored.cend(), input.begin());
            }
        }

        //- Set adaptor for different input, copying input as required
        void set(Container<InputType>& input, const bool doCopy = true)
        {
            if (orig_.get() != &input)
            {
                // Commit changes to old input first
                this->commit();
            }
            this->setInput(input, doCopy);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

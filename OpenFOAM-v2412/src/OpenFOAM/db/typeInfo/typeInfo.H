/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

Description
    Basic run-time type information using word as the type's name.
    Used to enhance the standard RTTI to cover I/O.

    The user can get the type's type name using the type info access function
    \code
        type()
    \endcode

    The isA functions:
    \code
        isA<Type>(obj)
        isA_constCast<Type>(obj)
    \endcode
    which return const or non-const pointers to the cast object,
    nullptr if cast is not possible (can be tested as a bool).

    The reference type cast template function:
    \code
        refCast<Type>(obj)
        refConstCast<Type>(obj)
    \endcode
    wraps dynamic_cast to handle failed casts and generate a FatalError.

\*---------------------------------------------------------------------------*/

#ifndef Foam_typeInfo_H
#define Foam_typeInfo_H

#include "error.H"
#include "className.H"
#include <typeinfo>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Declarations (for use in header files)

//- Declare a ClassNameNoDebug() with extra virtual type info
#define TypeNameNoDebug(TypeNameString)                                        \
    ClassNameNoDebug(TypeNameString);                                          \
    virtual const word& type() const { return typeName; }

//- Declare a ClassName() with extra virtual type info
#define TypeName(TypeNameString)                                               \
    ClassName(TypeNameString);                                                 \
    virtual const word& type() const { return typeName; }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Attempt dynamic_cast to \c Type.
//  \note The template types should \em not include any \c const qualifier.
//  \return const pointer to cast object, nullptr if cast is not possible
template<class Type, class U>
inline const Type* isA(const U& obj)
{
    const U* p = &obj;
    return dynamic_cast<const Type*>(p);
}


//- Attempt dynamic_cast to \c Type followed by a const_cast of the result.
//  \note The template types should \em not include any \c const qualifier.
//  \return non-const pointer to cast object, nullptr if cast is not possible
template<class Type, class U>
inline Type* isA_constCast(const U& obj)
{
    const U* p = &obj;
    return const_cast<Type*>(dynamic_cast<const Type*>(p));
}


//- Check if typeid of the object and \c Type are identical
template<class Type, class U>
inline bool isType(const U& obj)
{
    return typeid(Type) == typeid(obj);
}


//- A dynamic_cast (for references) that generates FatalError on failed casts
//  Respects the constness of the template types.
template<class Type, class U>
inline Type& dynamicCast(U& obj)
{
    U* p = &obj;
    Type* casted = dynamic_cast<Type*>(p);

    if (!casted)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << typeid(U).name()
            << " to type " << typeid(Type).name()
            << abort(FatalError);
    }

    return *casted;
}


//- A dynamic_cast (for references) that generates FatalIOError on failed casts
//  Respects the constness of the template types.
template<class Type, class U>
inline Type& dynamicCast(U& obj, const dictionary& dict)
{
    U* p = &obj;
    Type* casted = dynamic_cast<Type*>(p);

    if (!casted)
    {
        FatalIOErrorInFunction(dict)
            << "Attempt to cast type " << typeid(U).name()
            << " to type " << typeid(Type).name()
            << abort(FatalIOError);
    }

    return *casted;
}


//- A dynamic_cast (for references) to \c Type reference.
//  \note Respects the constness of the template types.
//  \return reference to cast object, or FatalError on failed casts
//      and use the virtual type() method for error messages.
template<class Type, class U>
inline Type& refCast(U& obj)
{
    U* p = &obj;
    Type* casted = dynamic_cast<Type*>(p);

    if (!casted)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << obj.type()
            << " to type " << Type::typeName
            << abort(FatalError);
    }

    return *casted;
}


//- A dynamic_cast (for const references) to \c Type reference,
//- followed by a const_cast of the result.
//  \note The template types should \em not include any \c const qualifier.
//  \return non-const reference to cast object, or FatalError on failed casts
//      and use the virtual type() method for error messages.
template<class Type, class U>
inline Type& refConstCast(const U& obj)
{
    const U* p = &obj;
    const Type* casted = dynamic_cast<const Type*>(p);

    if (!casted)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << obj.type()
            << " to type " << Type::typeName
            << abort(FatalError);
    }

    return const_cast<Type&>(*casted);
}


//- A dynamic_cast (for references) that generates FatalIOError on failed casts,
//- uses the virtual type() method for error messages.
//  Respects the constness of the template types.
template<class Type, class U>
inline Type& refCast(U& obj, const dictionary& dict)
{
    U* p = &obj;
    Type* casted = dynamic_cast<Type*>(p);

    if (!casted)
    {
        FatalIOErrorInFunction(dict)
            << "Attempt to cast type " << obj.type()
            << " to type " << Type::typeName
            << abort(FatalIOError);
    }

    return *casted;
}


//- A dynamic_cast (for references) that generates FatalError on failed casts,
//- uses the virtual type() method for error messages.
//- The index can be used to convey additional context.
//  Respects the constness of the template types.
template<class Type, class U>
inline Type& refCast(U& obj, const label index)
{
    U* p = &obj;
    Type* casted = dynamic_cast<Type*>(p);

    if (!casted)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << obj.type()
            << " to type " << Type::typeName
            << " at index " << index
            << abort(FatalError);
    }

    return *casted;
}


// * * * * * * * * * * * * * * * * Functors  * * * * * * * * * * * * * * * * //

//- Test if dynamic_cast to Type is possible, as a functor
template<class Type>
struct isAOp
{
    template<class U>
    bool operator()(const U& obj) const
    {
        return isA<Type, U>(obj);
    }
};


//- Test if typeid is identical to the Type, as a functor
template<class Type>
struct isTypeOp
{
    template<class U>
    bool operator()(const U& obj) const
    {
        return isType<Type, U>(obj);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

Namespace
    Foam::glTF

Description
    Namespace for handling glTF creation.
    https://www.khronos.org/registry/glTF/

Class
    Foam::glTF::base

Description
    Base class for glTF entities

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_base_H
#define foam_gltf_base_H

#include "word.H"
#include "label.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{

    enum class componentTypes : int
    {
        BYTE            = 5120, //!< 1 byte
        UNSIGNED_BYTE   = 5121, //!< 1 byte
        SHORT           = 5122, //!< 2 bytes
        UNSIGNED_SHORT  = 5123, //!< 2 bytes
        UNSIGNED_INT    = 5125, //!< 4 bytes
        FLOAT           = 5126  //!< 4 bytes
    };

    // accessor
    enum class dataTypes
    {
        SCALAR, //!< 1 component
        VEC2,   //!< 2 components
        VEC3,   //!< 3 components
        VEC4,   //!< 4 components
        MAT2,   //!< 4 components
        MAT3,   //!< 9 components
        MAT4    //!< 16 components
    };

    enum class attributeTypes
    {
        POSITION,   //!< VEC3 XYZ vertex positions; requires 'min' and 'max'
        NORMAL,     //!< VEC3 Normalised XYZ vertex normals
        TANGENT,    //!< VEC4 XYZW vertex tangents
        TEXCOORD_0, //!< VEC2 UV texture coordinates
        TEXCOORD_1, //!< VEC2 UV texture coordinates
        COLOR_0,    //!< VEC3 (rgb), VEC4 (rgba)
        JOINTS_0,   //!< VEC4
        WEIGHTS_0   //!< VEC4
    };

    // bufferView
    enum class targetTypes : int
    {
        ARRAY_BUFFER = 34962, //!< vertex attributes
        ELEMENT_ARRAY_BUFFER = 34963 //!< vertex indices
    };

    enum class primitiveModes : int
    {
        POINTS = 0,
        LINES = 1,
        LINE_LOOP = 2,
        LINE_STRIP = 3,
        TRIANGLES = 4,
        TRIANGLE_STRIP = 5,
        TRIANGLE_FAN = 6
    };

    // Helper function to return the enum value
    template<class Type>
    auto key(const Type& t) -> typename std::enable_if
    <
        std::is_enum<Type>::value,
        typename std::underlying_type<Type>::type
    >::type
    {
        return static_cast<typename std::underlying_type<Type>::type>(t);
    }

/*---------------------------------------------------------------------------*\
                          Class base Declaration
\*---------------------------------------------------------------------------*/

class base
{
protected:

    // Protected Data

        //- Name
        word name_;

        //- ID
        label id_;


public:

    // Constructors

        //- Default construct
        base()
        :
            id_(-1)
        {}

        //- Construct with name
        explicit base(const word& name)
        :
            name_(name),
            id_(-1)
        {}


    //- Destructor
    ~base() = default;


    // Member Functions

        //- Return access to the ID
        label& id() noexcept
        {
            return id_;
        }

        //- Return const access to the name
        const word& name() const noexcept
        {
            return name_;
        }

        //- Write
        void write(Ostream& os) const
        {
            os  << "," << nl
                << indent << "\"name\" : \"" << name_ << "\"";
        }


    // Member Operators

        void operator=(const base& gltf)
        {
            name_ = gltf.name_;
            id_ = gltf.id_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

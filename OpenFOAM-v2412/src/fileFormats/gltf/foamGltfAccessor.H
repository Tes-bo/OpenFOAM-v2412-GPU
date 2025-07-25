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

Class
    Foam::glTF::accessor

Description
    glTF accessor

Note
    Implements the glTF v2 specification

SourceFiles
    foamGltfAccessor.C
    foamGltfAccessorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_accessor_H
#define foam_gltf_accessor_H

#include "foamGltfBase.H"
#include "Field.H"
#include "label.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{
    class accessor;
}

Ostream& operator<<(Ostream& os, const glTF::accessor& acc);

namespace glTF
{

/*---------------------------------------------------------------------------*\
                       Class glTF::accessor Declaration
\*---------------------------------------------------------------------------*/

class accessor
:
    public base
{
protected:

    // Protected Data

        //- Buffer view index
        label bufferViewId_;

        //- Byte offset
        label byteOffset_;

        //- Component type
        label componentType_;

        //- Data size
        label count_;

        //- Data type
        string type_;

        //- Max value.  Note: stored as a string for convenience
        string max_;

        //- Min value.  Note: stored as a string for convenience
        string min_;

        //- Flag to indicate whether min and max values are available
        bool minMax_;


    // Protected Member Functions

        //- Return the glTF value type for the given OpenFOAM type
        template<class Type>
        static string getValueType();

        //- Stringify the value
        template<class Type>
        static string toString(const Type& val);


public:

    // Constructors

        //- Default construct
        accessor();

        //- Construct from name
        explicit accessor(const word& name);


    //- Destructor
    ~accessor() = default;


    // Public Member Functions

        //- Return the buffer view index
        label& bufferViewId() noexcept;

        //- Return the byte offset
        label& byteOffset() noexcept;

        //- Return the component type
        label& componentType() noexcept;

        //- Return the data size
        label& count() noexcept;

        //- Return the type
        string& type() noexcept;

        //- Set the accessor
        template<class Type>
        void set(const Field<Type>& fld, bool calcMinMax = true);

        //- Write
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamGltfAccessorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

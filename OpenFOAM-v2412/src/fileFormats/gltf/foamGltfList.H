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
    Foam::glTF::List

Description
    Container for glTF entities

Note
    Implements the glTF v2 specification

SourceFiles
    foamGltfList.C

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_list_H
#define foam_gltf_list_H

#include "DynamicList.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{
    template<class Type> class List;
}

template<class Type>
Ostream& operator<<(Ostream& os, const glTF::List<Type>& lst);

namespace glTF
{

/*---------------------------------------------------------------------------*\
                         Class glTF::List Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class List
{
    // Private Data

        //- Storage
        DynamicList<Type> data_;


public:

    // Constructors

        //- Default construct
        List() = default;


    //- Destructor
    ~List() = default;


    // Public Member Functions

        //- Helper to create a new Type on the list and set the ID
        Type& create(const word& name);

        //- Return const access to the underlying list
        const DynamicList<Type>& data() const noexcept;

        //- List contains no data
        bool empty() const noexcept;

        //- Return the list size
        label size() const noexcept;

        //- Write
        void write(Ostream& os, const word& keyword, bool firstEntry = false);

        //- Write
        void write(Ostream& os) const;

        Type& operator[](const label i);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamGltfList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

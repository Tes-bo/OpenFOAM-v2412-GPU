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
    Foam::glTF::object

Description
    glTF binary object

Note
    Implements the glTF v2 specification

SourceFiles
    foamGltfObject.C
    foamGltfObjectTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_object_H
#define foam_gltf_object_H

#include "foamGltfBase.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{

/*---------------------------------------------------------------------------*\
                        Class glTF::object Declaration
\*---------------------------------------------------------------------------*/

class object
:
    public base
{
    // Private Data

        //- Buffer storage
        Foam::List<float> data_;


public:

    // Constructors

        //- Default construct
        object();

        //- Construct with name
        explicit object(const word& name);


    // Public Member Functions

        //- Add data to the buffer
        template<class Type>
        void addData(const Type& fld);

        //- Add data to the buffer from 2 containers of the same size
        //  E.g. to combine vector and scalar to create RGBA data
        template<class Type1, class Type2>
        void addData(const Type1& fld1, const Type2& fld2);

        //- Return const access to the data buffer
        const Foam::List<float>& data() const noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamGltfObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

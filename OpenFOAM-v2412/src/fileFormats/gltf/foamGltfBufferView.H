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
    Foam::glTF::bufferView

Description
    glTF buffer view - provides a view/slice of the glTF buffer

Note
    Implements the glTF v2 specification

SourceFiles
    foamGltfBufferView.C

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_bufferview_H
#define foam_gltf_bufferview_H

#include "foamGltfBase.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{
    class bufferView;
}

Ostream& operator<<(Ostream&, const glTF::bufferView& bv);

namespace glTF
{

/*---------------------------------------------------------------------------*\
                      Class glTF::bufferView Declaration
\*---------------------------------------------------------------------------*/

class bufferView
:
    public base
{
protected:

    // Protected Data

        //- Buffer index
        label buffer_;

        //- Byte offset
        label byteOffset_;

        //- Byte length
        label byteLength_;

        //- Target
        label target_;


public:

    // Constructors

        //- Default construct
        bufferView();

        //- Construct with name
        explicit bufferView(const word& name);


    //- Destructor
    ~bufferView() = default;


    // Public Member Functions

        //- Return the buffer index
        label& buffer() noexcept;

        //- Return the byte offset
        label& byteOffset() noexcept;

        //- Return the byte length
        label& byteLength() noexcept;

        //- Return the target
        label& target() noexcept;

        //- Write
        void write(Ostream& os) const;


    // Member Operators

        void operator=(const bufferView& bv);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

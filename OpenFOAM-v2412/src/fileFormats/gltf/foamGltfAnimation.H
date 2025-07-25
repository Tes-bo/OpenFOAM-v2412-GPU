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
    Foam::glTF::animation

Description
    glTF animation

Note
    Implements the glTF v2 specification

SourceFiles
    foamGltfAnimation.C

\*---------------------------------------------------------------------------*/

#ifndef foam_gltf_animation_H
#define foam_gltf_animation_H

#include "foamGltfBase.H"
#include "DynamicList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace glTF
{
    class animation;
}

Ostream& operator<<(Ostream& os, const glTF::animation& animation);

namespace glTF
{

/*---------------------------------------------------------------------------*\
                       Class glTF::animation Declaration
\*---------------------------------------------------------------------------*/


class animation
:
    public base
{
protected:

    // Local Helpers

        // Sampler
        struct glTFSampler
        {
            label input;
            string interpolation;
            label output;
        };

        // Channel target
        struct glTFTarget
        {
            label node;
            string path;
        };

        // Channel
        struct glTFChannel
        {
            label samplerId;
            glTFTarget target;
        };



    // Protected Member Data

        //- Samplers
        DynamicList<glTFSampler> samplers_;

        //- Channels
        DynamicList<glTFChannel> channels_;


public:

    // Constructors

        //- Default construct
        animation();

        //- Construct with name
        explicit animation(const word& name);


    //- Destructor
    ~animation() = default;


    // Public Member Functions

        //- Add translation
        void addTranslation
        (
            const label inputId,
            const label outputId,
            const label nodeId,
            const string& interpolation
        );

        //- Write
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace glTF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

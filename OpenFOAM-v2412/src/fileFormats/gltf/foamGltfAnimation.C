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

\*---------------------------------------------------------------------------*/

#include "foamGltfAnimation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::glTF::animation::animation()
:
    base(),
    samplers_(),
    channels_()
{}


Foam::glTF::animation::animation(const word& name)
:
    base(name),
    samplers_(),
    channels_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::glTF::animation::addTranslation
(
    const label inputId,
    const label outputId,
    const label nodeId,
    const string& interpolation
)
{
    glTFSampler sampler;
    sampler.input = inputId;
    sampler.output = outputId;
    sampler.interpolation = interpolation;
    samplers_.append(sampler);

    glTFChannel channel;
    channel.samplerId = samplers_.size() - 1;
    channel.target.node = nodeId;
    channel.target.path = "translation";
    channels_.append(channel);
}


void Foam::glTF::animation::write(Ostream& os) const
{
    os  << indent << "\"samplers\" : [" << nl << incrIndent;

    forAll(samplers_, i)
    {
        const auto& sampler = samplers_[i];

        os  << indent << "{" << nl << incrIndent
            << indent << "\"input\" : " << sampler.input << "," << nl
            << indent << "\"interpolation\" : " << sampler.interpolation
            << "," << nl
            << indent << "\"output\" : " << sampler.output << nl
            << decrIndent << indent << "}";

        if (i != samplers_.size() - 1) os  << "," << nl;
    }

    os  << nl << decrIndent << indent << "]," << nl;

    os  << indent << "\"channels\" : [" << nl << incrIndent;

    forAll(channels_, i)
    {
        const auto& channel = channels_[i];

        os  << indent << "{" << nl << incrIndent
            << indent << "\"sampler\" : " << channel.samplerId << "," << nl
            << indent << "\"target\" : {" << incrIndent << nl
            << indent << "\"node\" : " << channel.target.node << "," << nl
            << indent << "\"path\" : " << channel.target.path << nl
            << decrIndent << indent << "}" << nl
            << decrIndent << indent << "}";

        if (i != channels_.size() - 1) os  << "," << nl;
    }

    os  << nl << decrIndent << indent << "]";
}


Foam::Ostream& Foam::operator<<(Ostream& os, const glTF::animation& animation)
{
    animation.write(os);

    return os;
}


// ************************************************************************* //

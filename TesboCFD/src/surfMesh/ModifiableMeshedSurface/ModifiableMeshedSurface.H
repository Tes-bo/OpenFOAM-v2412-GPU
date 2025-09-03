/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::ModifiableMeshedSurface

Description
    A special purpose MeshedSurface that exposes the stored values
    for direct modification.

    Its usage should be restricted to special cases where the surface
    needs modifications as an atomic operation.

\*---------------------------------------------------------------------------*/

#ifndef ModifiableMeshedSurface_H
#define ModifiableMeshedSurface_H

#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class ModifiableMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class ModifiableMeshedSurface
:
    public MeshedSurface<Face>
{
public:

    // Constructors

        //- Default construct. Use swap or transfer to manage content
        ModifiableMeshedSurface()
        :
            MeshedSurface<Face>()
        {}


        //- Construct by transferring the contents from a MeshedSurface
        explicit ModifiableMeshedSurface(MeshedSurface<Face>&& surf)
        :
            MeshedSurface<Face>(std::move(surf))
        {}


    //- Destructor
    virtual ~ModifiableMeshedSurface() = default;


    // Member Functions

        // Expose various protected methods
        using MeshedSurface<Face>::storedFaces;
        using MeshedSurface<Face>::storedPoints;
        using MeshedSurface<Face>::storedZones;

        //- Copy and Move assignments
        using MeshedSurface<Face>::operator=;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

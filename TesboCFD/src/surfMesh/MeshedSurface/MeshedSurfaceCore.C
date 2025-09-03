/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "MeshedSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{

    // Transcribe 'face' to 'face' (ie, just transfer)
    template<>
    void Foam::MeshedSurface<Foam::face>::transcribe
    (
        MeshedSurface<face>& surf
    )
    {
        this->transfer(surf);
        this->addZonesToFaces(); // currently a no-op
    }

    // Transcribe 'face' to 'triFace'
    // Transfer points/zones and triangulate faces
    template<>
    void Foam::MeshedSurface<Foam::triFace>::transcribe
    (
        MeshedSurface<face>& surf
    )
    {
        // First triangulate.
        // Potentially wasteful of space, but adjusts zones and
        // invalidates the faceIds
        surf.triangulate();

        this->storedPoints().transfer(surf.storedPoints());
        this->storedZones().transfer(surf.storedZones());

        // Transcribe from face -> triFace
        const List<face>& origFaces = surf.surfFaces();
        List<triFace> newFaces(origFaces.size());
        forAll(origFaces, facei)
        {
            newFaces[facei] = triFace
            (
                static_cast<const labelUList&>(origFaces[facei])
            );
        }
        surf.clear();

        this->storedFaces().transfer(newFaces);
        this->addZonesToFaces(); // currently a no-op
    }


    // Transcribe 'face' to 'labelledTri'
    // Transfer points/zones and triangulate faces
    template<>
    void Foam::MeshedSurface<Foam::labelledTri>::transcribe
    (
        MeshedSurface<face>& surf
    )
    {
        // First triangulate.
        // Potentially wasteful of space, but adjusts zones and
        // invalidates the faceIds
        surf.triangulate();

        this->storedPoints().transfer(surf.storedPoints());
        this->storedZones().transfer(surf.storedZones());

        // Transcribe from face -> labelledTri (via triFace)
        const List<face>& origFaces = surf.surfFaces();
        List<labelledTri> newFaces(origFaces.size());
        forAll(origFaces, facei)
        {
            newFaces[facei] = triFace
            (
                static_cast<const labelUList&>(origFaces[facei])
            );
        }
        surf.clear();

        this->storedFaces().transfer(newFaces);
        this->addZonesToFaces(); // for labelledTri
    }


    // Propagate zone information on face regions for labelledTri.
    template<>
    bool Foam::MeshedSurface<Foam::labelledTri>::addZonesToFaces()
    {
        List<labelledTri>& faceLst = this->storedFaces();
        const surfZoneList& zones = this->surfZones();

        forAll(zones, zoneI)
        {
            const surfZone& zone = zones[zoneI];

            label faceI = zone.start();
            forAll(zone, i)
            {
                faceLst[faceI++].region() = zoneI;
            }
        }

        return true;
    }


}  // End namespace Foam


// ************************************************************************* //

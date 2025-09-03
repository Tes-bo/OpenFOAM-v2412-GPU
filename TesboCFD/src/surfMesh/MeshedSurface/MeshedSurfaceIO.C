/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
#include "boundBox.H"
#include "faceTraits.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
Foam::Istream& Foam::MeshedSurface<Face>::read(Istream& is)
{
    is  >> this->storedZones()
        >> this->storedPoints()
        >> this->storedFaces();

    is.check(FUNCTION_NAME);
    return is;
}


template<class Face>
Foam::Ostream& Foam::MeshedSurface<Face>::write(Ostream& os) const
{
    os  << this->surfZones()
        << this->points()
        << this->surfFaces();

    os.check(FUNCTION_NAME);
    return os;
}


template<class Face>
void Foam::MeshedSurface<Face>::writeStats(Ostream& os) const
{
    os  << "points      : " << this->points().size() << nl;
    if (faceTraits<Face>::isTri())
    {
        os << "triangles   : " << this->size() << nl;
    }
    else
    {
        label nTri = 0, nQuad = 0;
        for (const Face& f : *this)
        {
            const label n = f.size();

            if (n == 3)
            {
                ++nTri;
            }
            else if (n == 4)
            {
                ++nQuad;
            }
        }

        os  << "faces       : " << this->size()
            << "  (tri:" << nTri << " quad:" << nQuad
            << " poly:" << (this->size() - nTri - nQuad) << ")" << nl;
    }

    os  << "boundingBox : " << boundBox(this->points()) << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Face>
Foam::Istream& Foam::operator>>(Istream& is, MeshedSurface<Face>& surf)
{
    return surf.read(is);
}


template<class Face>
Foam::Ostream& Foam::operator<<(Ostream& os, const MeshedSurface<Face>& surf)
{
    return surf.write(os);
}


// ************************************************************************* //

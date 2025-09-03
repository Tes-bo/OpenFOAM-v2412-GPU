/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

template<class FaceType>
void Foam::meshTools::writeOBJ
(
    Ostream& os,
    const UList<FaceType>& faces,
    const UList<point>& points,
    const labelList& faceLabels
)
{
    Map<label> foamToObj(4*faceLabels.size());

    label vertI = 0;

    for (const label facei : faceLabels)
    {
        const FaceType& f = faces[facei];

        forAll(f, fp)
        {
            if (foamToObj.insert(f[fp], vertI))
            {
                writeOBJ(os, points[f[fp]]);
                vertI++;
            }
        }

        os << 'f';
        forAll(f, fp)
        {
            os << ' ' << foamToObj[f[fp]]+1;
        }
        os << ' ' << foamToObj[f[0]]+1 << endl;
    }
}


template<class FaceType>
void Foam::meshTools::writeOBJ
(
    Ostream& os,
    const UList<FaceType>& faces,
    const UList<point>& points
)
{
    labelList allFaces(faces.size());
    forAll(allFaces, i)
    {
        allFaces[i] = i;
    }
    writeOBJ(os, faces, points, allFaces);
}


// ************************************************************************* //

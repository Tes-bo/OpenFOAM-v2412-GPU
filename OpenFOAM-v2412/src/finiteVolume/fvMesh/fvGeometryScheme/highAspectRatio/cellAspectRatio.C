/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "cellAspectRatio.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellAspectRatio, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::cellAspectRatio::cellAspectRatio(const polyMesh& mesh)
:
    MeshObject_type(mesh)
{
    calcAspectRatio();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::cellAspectRatio::~cellAspectRatio()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellAspectRatio::calcAspectRatio()
{
    if (debug)
    {
        InfoInFunction << "Calculating cell aspect ratio" << endl;
    }

    const polyMesh& mesh = mesh_;
    const pointField& cellCentres = mesh.cellCentres();
    const scalarField& cellVolumes = mesh.cellVolumes();
    const vectorField& faceAreas = mesh.faceAreas();
    const vectorField& faceCentres = mesh.faceCentres();
    const cellList& cells = mesh.cells();
    //const faceList& faces = mesh.faces();
    //const pointField& points = mesh.points();

    scalarField& aRatio = *this;
    aRatio.setSize(mesh.nCells());

    forAll(cells, celli)
    {
        const point& cc = cellCentres[celli];
        const cell& cFaces = cells[celli];

        scalar sumA = Zero;
        scalar maxMag = Zero;

        for (const label facei : cFaces)
        {
            const vector& n = faceAreas[facei];

            sumA += mag(n);

            //// Max distance from point to cell centre
            //const face& f = faces[facei];
            //for (const label pointi : f)
            //{
            //    const point& pt = points[pointi];
            //    const vector d(pt-cc);
            //    maxMag = max(maxMag, magSqr(d));
            //}

            // Max distance from face centre to cell centre
            const point& fc = faceCentres[facei];
            maxMag = max(maxMag, magSqr(fc-cc));
        }
        sumA /= cFaces.size();

        aRatio[celli] = 1.0;
        if (sumA > ROOTVSMALL)
        {
            // Local length scale
            const scalar length = cellVolumes[celli]/sumA;

            if (length > ROOTVSMALL)
            {
                // Max edge length
                maxMag = Foam::sqrt(maxMag);

                //aRatio[celli] = Foam::sqrt(4.0/3.0)*maxMag/length;
                aRatio[celli] = 2.0*maxMag/length;
            }
        }
    }

    if (debug)
    {
        InfoInFunction << "Calculated cell aspect ratio min:" << gMin(aRatio)
            << " max:" << gMax(aRatio) << " average:" << gAverage(aRatio)
            << endl;
    }
}


// ************************************************************************* //

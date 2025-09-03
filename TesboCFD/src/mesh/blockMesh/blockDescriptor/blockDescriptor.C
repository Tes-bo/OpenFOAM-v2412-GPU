/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "blockDescriptor.H"
#include "blockMeshTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::blockDescriptor::assignGradings
(
    const UList<gradingDescriptors>& ratios
)
{
    bool ok = true;

    switch (ratios.size())
    {
        case 0:
        {
            expand_.resize(12);
            expand_ = gradingDescriptors();
            break;
        }
        case 1:
        {
            // Identical in x/y/z-directions
            expand_.resize(12);
            expand_ = ratios[0];
            break;
        }
        case 3:
        {
            expand_.resize(12);

            // x-direction
            expand_[0]  = ratios[0];
            expand_[1]  = ratios[0];
            expand_[2]  = ratios[0];
            expand_[3]  = ratios[0];

            // y-direction
            expand_[4]  = ratios[1];
            expand_[5]  = ratios[1];
            expand_[6]  = ratios[1];
            expand_[7]  = ratios[1];

            // z-direction
            expand_[8]  = ratios[2];
            expand_[9]  = ratios[2];
            expand_[10] = ratios[2];
            expand_[11] = ratios[2];
            break;
        }
        case 12:
        {
            expand_ = ratios;
            break;
        }
        default:
        {
            ok = false;
            break;
        }
    }

    return ok;
}


void Foam::blockDescriptor::check(const Istream& is)
{
    for (const label pointi : blockShape_)
    {
        if (pointi < 0 || pointi >= vertices_.size())
        {
            FatalIOErrorInFunction(is)
                << "Point label (" << pointi
                << ") out of range 0.." << vertices_.size() - 1
                << " in block " << *this
                << exit(FatalIOError);
        }
    }

    const point blockCentre(blockShape_.centre(vertices_));
    const faceList faces(blockShape_.faces());

    // Check each face is outward-pointing with respect to the block centre
    label outwardFaceCount = 0;
    boolList correctFaces(faces.size(), true);

    forAll(faces, i)
    {
        point faceCentre(faces[i].centre(vertices_));
        vector faceNormal(faces[i].areaNormal(vertices_));

        if (mag(faceNormal) > SMALL)
        {
            if (((faceCentre - blockCentre) & faceNormal) > 0)
            {
                outwardFaceCount++;
            }
            else
            {
                correctFaces[i] = false;
            }
        }
        else
        {
            outwardFaceCount++;
        }
    }

    // If all faces are inward-pointing the block is inside-out
    if (outwardFaceCount == 0)
    {
        FatalIOErrorInFunction(is)
            << "Block " << *this << " is inside-out"
            << exit(FatalIOError);
    }
    else if (outwardFaceCount != faces.size())
    {
        FatalIOErrorInFunction(is)
            << "Block " << *this << " has inward-pointing faces"
            << nl << "    ";

        forAll(correctFaces, i)
        {
            if (!correctFaces[i])
            {
                FatalIOError<< faces[i] << token::SPACE;
            }
        }

        FatalIOError << exit(FatalIOError);
    }
}


void Foam::blockDescriptor::findCurvedFaces(const label blockIndex)
{
    const faceList shapeFaces(blockShape().faces());

    forAll(shapeFaces, shapeFacei)
    {
        forAll(blockFaces_, facei)
        {
            const face& f = blockFaces_[facei].vertices();

            // Accept (<block> <face>) face description
            if
            (
                (
                    f.size() == 2
                 && f[0] == blockIndex
                 && f[1] == shapeFacei
                )
             || face::sameVertices(f, shapeFaces[shapeFacei])
            )
            {
                curvedFaces_[shapeFacei] = facei;
                ++nCurvedFaces_;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockDescriptor::blockDescriptor
(
    const cellShape& bshape,
    const pointField& vertices,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    const labelVector& density,
    const UList<gradingDescriptors>& expand,
    const word& zoneName
)
:
    ijkMesh(density),
    vertices_(vertices),
    blockEdges_(edges),
    blockFaces_(faces),
    blockShape_(bshape),
    expand_(),
    index_(-1),
    zoneName_(zoneName),
    curvedFaces_(-1),
    nCurvedFaces_(0)
{
    if (!assignGradings(expand))
    {
        FatalErrorInFunction
            << "Unknown definition of expansion ratios: " << expand
            << exit(FatalError);
    }

    findCurvedFaces();
}


Foam::blockDescriptor::blockDescriptor
(
    const dictionary& dict,
    const label blockIndex,
    const pointField& vertices,
    const blockEdgeList& edges,
    const blockFaceList& faces,
    Istream& is
)
:
    ijkMesh(),
    vertices_(vertices),
    blockEdges_(edges),
    blockFaces_(faces),
    blockShape_(),
    expand_(),
    index_(blockIndex),
    zoneName_(),
    curvedFaces_(-1),
    nCurvedFaces_(0)
{
    // Read cell model and list of vertices (potentially with variables)
    word model(is);
    blockShape_ = cellShape
    (
        model,
        blockMeshTools::read<label>
        (
            is,
            dict.subOrEmptyDict("namedVertices")
        )
    );

    // Examine next token
    token t(is);

    // Optional zone name
    if (t.isWord())
    {
        zoneName_ = t.wordToken();

        // Examine next token
        is >> t;
    }
    is.putBack(t);

    if (t.isPunctuation())
    {
        // New-style: read a list of 3 values
        if (t.pToken() == token::BEGIN_LIST)
        {
            is >> ijkMesh::sizes();
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Incorrect token while reading n, expected '(', found "
                << t.info()
                << exit(FatalIOError);
        }
    }
    else
    {
        // Old-style: read three labels
        IOWarningInFunction(is)
            << "Encountered old-style specification of mesh divisions"
            << endl;

        is  >> ijkMesh::sizes().x()
            >> ijkMesh::sizes().y()
            >> ijkMesh::sizes().z();
    }

    is >> t;
    if (!t.isWord())
    {
        is.putBack(t);
    }

    List<gradingDescriptors> expand(is);

    if (!assignGradings(expand))
    {
        FatalErrorInFunction
            << "Unknown definition of expansion ratios: " << expand
            << exit(FatalError);
    }

    check(is);

    findCurvedFaces(blockIndex);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::FixedList<Foam::pointField, 6>
Foam::blockDescriptor::facePoints(const pointField& points) const
{
    const label ni = sizes().x();
    const label nj = sizes().y();
    const label nk = sizes().z();

    // Caches points for curvature correction
    FixedList<pointField, 6> facePoints;

    facePoints[0].setSize((nj + 1)*(nk + 1));
    facePoints[1].setSize((nj + 1)*(nk + 1));

    for (label j=0; j<=nj; j++)
    {
        for (label k=0; k<=nk; k++)
        {
            facePoints[0][facePointLabel(0, j, k)] =
                points[pointLabel(0, j, k)];
            facePoints[1][facePointLabel(1, j, k)] =
                points[pointLabel(ni, j, k)];
        }
    }

    facePoints[2].setSize((ni + 1)*(nk + 1));
    facePoints[3].setSize((ni + 1)*(nk + 1));

    for (label i=0; i<=ni; i++)
    {
        for (label k=0; k<=nk; k++)
        {
            facePoints[2][facePointLabel(2, i, k)] =
                points[pointLabel(i, 0, k)];
            facePoints[3][facePointLabel(3, i, k)] =
                points[pointLabel(i, nj, k)];
        }
    }

    facePoints[4].setSize((ni + 1)*(nj + 1));
    facePoints[5].setSize((ni + 1)*(nj + 1));

    for (label i=0; i<=ni; i++)
    {
        for (label j=0; j<=nj; j++)
        {
            facePoints[4][facePointLabel(4, i, j)] =
                points[pointLabel(i, j, 0)];
            facePoints[5][facePointLabel(5, i, j)] =
                points[pointLabel(i, j, nk)];
        }
    }

    return facePoints;
}


void Foam::blockDescriptor::correctFacePoints
(
    FixedList<pointField, 6>& facePoints
) const
{
    forAll(curvedFaces_, blockFacei)
    {
        if (curvedFaces_[blockFacei] >= 0)
        {
            blockFaces_[curvedFaces_[blockFacei]].project
            (
                *this,
                blockFacei,
                facePoints[blockFacei]
            );
        }
    }
}


void Foam::blockDescriptor::write
(
    Ostream& os,
    const label val,
    const dictionary& d
)
{
    const dictionary* varDictPtr = d.findDict("namedBlocks");
    if (varDictPtr)
    {
        blockMeshTools::write(os, val, *varDictPtr);
    }
    else
    {
        os << val;
    }
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const blockDescriptor& bd)
{
    const cellShape& bshape = bd.blockShape();
    const labelList& blockLabels = bshape;

    os  << bshape.model().name() << " (";

    forAll(blockLabels, labeli)
    {
        if (labeli)
        {
            os  << ' ';
        }
        os  << blockLabels[labeli];
    }
    os  << ')';

    if (bd.zoneName().size())
    {
        os  << ' ' << bd.zoneName();
    }

    os  << ' '  << bd.density()
        << " grading (";


    const List<gradingDescriptors>& expand = bd.grading();

    // Can we use a compact notation?
    if
    (
        // x-direction
        (
            expand[0] == expand[1]
         && expand[0] == expand[2]
         && expand[0] == expand[3]
        )
     && // y-direction
        (
            expand[4] == expand[5]
         && expand[4] == expand[6]
         && expand[4] == expand[7]
        )
     && // z-direction
        (
            expand[8] == expand[9]
         && expand[8] == expand[10]
         && expand[8] == expand[11]
        )
    )
    {
        os  << expand[0] << ' ' << expand[4] << ' ' << expand[8];
    }
    else
    {
        forAll(expand, edgei)
        {
            if (edgei)
            {
                os  << ' ';
            }
            os  << expand[edgei];
        }
    }

    os  << ')';

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<blockDescriptor>& iproxy
)
{
    const auto& bd = *iproxy;

    os  << "Dimensions:" << bd.density()
        << " nPoints:" << bd.nPoints()
        << " nCells:" << bd.nCells()
        << " nFaces:" << bd.nFaces()
        << " nInternalFaces:" << bd.nInternalFaces()
        << nl;

    return os;
}


// ************************************************************************* //

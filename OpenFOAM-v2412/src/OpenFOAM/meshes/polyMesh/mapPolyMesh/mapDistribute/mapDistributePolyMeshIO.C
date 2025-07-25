/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "mapDistributePolyMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistributePolyMesh::mapDistributePolyMesh
(
    const dictionary& dict,
    const label comm
)
:
    mapDistributePolyMesh(comm)
{
    mapDistributePolyMesh::readDict(dict);
}


Foam::mapDistributePolyMesh::mapDistributePolyMesh(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mapDistributePolyMesh::readDict(const dictionary& dict)
{
    // Cell information
    {
        const dictionary& subdict = dict.subDict("cellMap");

        subdict.readEntry("oldSize", nOldCells_);
        cellMap_.readDict(subdict);
    }

    // Face information
    {
        const dictionary& subdict = dict.subDict("faceMap");

        subdict.readEntry("oldSize", nOldFaces_);
        faceMap_.readDict(subdict);
    }

    // Point information
    {
        const dictionary& subdict = dict.subDict("pointMap");

        subdict.readEntry("oldSize", nOldPoints_);
        pointMap_.readDict(subdict);
    }

    // Patch information
    {
        const dictionary& subdict = dict.subDict("patchMap");

        subdict.readEntry("oldSizes", oldPatchSizes_);
        subdict.readEntry("oldStarts", oldPatchStarts_);
        subdict.readEntry("oldPointSizes", oldPatchNMeshPoints_);
        patchMap_.readDict(subdict);
    }
}


void Foam::mapDistributePolyMesh::writeCellMapEntries(Ostream& os) const
{
    os.beginBlock("cellMap");
    os.writeEntry("oldSize", nOldCells_);
    cellMap_.writeEntries(os);
    os.endBlock();
}


void Foam::mapDistributePolyMesh::writeFaceMapEntries(Ostream& os) const
{
    os.beginBlock("faceMap");
    os.writeEntry("oldSize", nOldFaces_);
    faceMap_.writeEntries(os);
    os.endBlock();
}


void Foam::mapDistributePolyMesh::writePointMapEntries(Ostream& os) const
{
    os.beginBlock("pointMap");
    os.writeEntry("oldSize", nOldPoints_);
    pointMap_.writeEntries(os);
    os.endBlock();
}


void Foam::mapDistributePolyMesh::writePatchMapEntries(Ostream& os) const
{
    os.beginBlock("patchMap");
    oldPatchSizes_.writeEntry("oldSizes", os);
    oldPatchStarts_.writeEntry("oldStarts", os);
    oldPatchNMeshPoints_.writeEntry("oldPointSizes", os);
    patchMap_.writeEntries(os);
    os.endBlock();
}


void Foam::mapDistributePolyMesh::writeEntries(Ostream& os) const
{
    writeCellMapEntries(os);

    os << nl;
    writeFaceMapEntries(os);

    os << nl;
    writePointMapEntries(os);

    os << nl;
    writePatchMapEntries(os);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, mapDistributePolyMesh& map)
{
    is.fatalCheck(FUNCTION_NAME);

    is  >> map.nOldPoints_
        >> map.nOldFaces_
        >> map.nOldCells_

        >> map.oldPatchSizes_
        >> map.oldPatchStarts_
        >> map.oldPatchNMeshPoints_

        >> map.pointMap_
        >> map.faceMap_
        >> map.cellMap_
        >> map.patchMap_;

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const mapDistributePolyMesh& map)
{
    os  << map.nOldPoints_ << token::SPACE
        << map.nOldFaces_ << token::SPACE
        << map.nOldCells_ << token::NL

        << map.oldPatchSizes_ << token::NL
        << map.oldPatchStarts_ << token::NL
        << map.oldPatchNMeshPoints_ << token::NL

        << map.pointMap_ << token::NL
        << map.faceMap_ << token::NL
        << map.cellMap_ << token::NL
        << map.patchMap_;

    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<mapDistributePolyMesh>& iproxy
)
{
    const auto& map = *iproxy;

    os.beginBlock("cellMap");
    os.writeEntry("oldSize", map.nOldCells());
    os  << map.cellMap().info();
    os.endBlock();

    os.beginBlock("faceMap");
    os.writeEntry("oldSize", map.nOldFaces());
    os  << map.faceMap().info();
    os.endBlock();

    os.beginBlock("pointMap");
    os.writeEntry("oldSize", map.nOldPoints());
    os  << map.pointMap().info();
    os.endBlock();

    return os;
}


// ************************************************************************* //

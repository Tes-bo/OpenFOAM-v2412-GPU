/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "IOobject.H"
#include "List.H"
#include "hexRef8Data.H"
#include "mapPolyMesh.H"
#include "mapDistributePolyMesh.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "refinementHistory.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hexRef8Data::hexRef8Data(const IOobject& io)
{
    if
    (
        // Or: !io.anyRead() ?
        io.readOpt() != IOobject::MUST_READ
     && io.readOpt() != IOobject::LAZY_READ
    )
    {
        return;
    }

    {
        typedef labelIOList Type;
        IOobject rio(io, "cellLevel");

        // haveFile
        if (returnReduceOr(rio.typeHeaderOk<Type>(true)))
        {
            Info<< "Reading hexRef8 data : " << rio.name() << endl;
            cellLevelPtr_.reset(new Type(rio));
        }
    }
    {
        typedef labelIOList Type;
        IOobject rio(io, "pointLevel");

        // haveFile
        if (returnReduceOr(rio.typeHeaderOk<Type>(true)))
        {
            Info<< "Reading hexRef8 data : " << rio.name() << endl;
            pointLevelPtr_.reset(new Type(rio));
        }
    }
    {
        typedef uniformDimensionedScalarField Type;
        IOobject rio(io, "level0Edge");

        // haveFile
        if (returnReduceOr(rio.typeHeaderOk<Type>(true)))
        {
            Info<< "Reading hexRef8 data : " << rio.name() << endl;
            level0EdgePtr_.reset(new Type(rio));
        }
    }
    {
        typedef refinementHistory Type;
        IOobject rio(io, "refinementHistory");

        // haveFile
        if (returnReduceOr(rio.typeHeaderOk<Type>(true)))
        {
            Info<< "Reading hexRef8 data : " << rio.name() << endl;
            refHistoryPtr_.reset(new Type(rio));
        }
    }
}


Foam::hexRef8Data::hexRef8Data
(
    const IOobject& io,
    const hexRef8Data& data,
    const labelList& cellMap,
    const labelList& pointMap
)
{
    if (data.cellLevelPtr_)
    {
        IOobject rio(io, data.cellLevelPtr_().name());

        cellLevelPtr_.reset
        (
            new labelIOList
            (
                rio,
                labelUIndList(data.cellLevelPtr_(), cellMap)()
            )
        );
    }
    if (data.pointLevelPtr_)
    {
        IOobject rio(io, data.pointLevelPtr_().name());

        pointLevelPtr_.reset
        (
            new labelIOList
            (
                rio,
                labelUIndList(data.pointLevelPtr_(), pointMap)()
            )
        );
    }
    if (data.level0EdgePtr_)
    {
        IOobject rio(io, data.level0EdgePtr_().name());

        level0EdgePtr_.reset
        (
            new uniformDimensionedScalarField(rio, data.level0EdgePtr_())
        );
    }
    if (data.refHistoryPtr_)
    {
        IOobject rio(io, data.refHistoryPtr_().name());

        refHistoryPtr_ = data.refHistoryPtr_().clone(rio, cellMap);
    }
}


Foam::hexRef8Data::hexRef8Data
(
    const IOobject& io,
    const UPtrList<const labelList>& cellMaps,
    const UPtrList<const labelList>& pointMaps,
    const UPtrList<const hexRef8Data>& procDatas
)
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(io.db());

    // cellLevel

    if (procDatas[0].cellLevelPtr_)
    {
        IOobject rio(io, procDatas[0].cellLevelPtr_().name());

        cellLevelPtr_.reset(new labelIOList(rio, mesh.nCells()));
        auto& cellLevel = *cellLevelPtr_;

        forAll(procDatas, procI)
        {
            const labelList& procCellLevel = procDatas[procI].cellLevelPtr_();
            labelUIndList(cellLevel, cellMaps[procI]) = procCellLevel;
        }
    }


    // pointLevel

    if (procDatas[0].pointLevelPtr_)
    {
        IOobject rio(io, procDatas[0].pointLevelPtr_().name());

        pointLevelPtr_.reset(new labelIOList(rio, mesh.nPoints()));
        auto& pointLevel = *pointLevelPtr_;

        forAll(procDatas, procI)
        {
            const labelList& procPointLevel = procDatas[procI].pointLevelPtr_();
            labelUIndList(pointLevel, pointMaps[procI]) = procPointLevel;
        }
    }


    // level0Edge

    if (procDatas[0].level0EdgePtr_)
    {
        IOobject rio(io, procDatas[0].level0EdgePtr_().name());

        level0EdgePtr_.reset
        (
            new uniformDimensionedScalarField
            (
                rio,
                procDatas[0].level0EdgePtr_()
            )
        );
    }


    // refinementHistory

    if (procDatas[0].refHistoryPtr_)
    {
        IOobject rio(io, procDatas[0].refHistoryPtr_().name());

        UPtrList<const refinementHistory> procRefs(procDatas.size());
        forAll(procDatas, i)
        {
            procRefs.set(i, &procDatas[i].refHistoryPtr_());
        }

        refHistoryPtr_.reset
        (
            new refinementHistory
            (
                rio,
                cellMaps,
                procRefs
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hexRef8Data::~hexRef8Data()
{}  // refinementHistory forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hexRef8Data::sync(const IOobject& io)
{
    const polyMesh& mesh = dynamic_cast<const polyMesh&>(io.db());

    if (returnReduceOr(cellLevelPtr_) && !cellLevelPtr_)
    {
        IOobject rio(io, "cellLevel");
        rio.readOpt(IOobject::NO_READ);
        cellLevelPtr_.reset
        (
            new labelIOList(rio, labelList(mesh.nCells(), Zero))
        );
    }
    // Take over IOobject settings
    if (cellLevelPtr_)
    {
        cellLevelPtr_->IOobjectOption::operator=(io);
    }

    if (returnReduceOr(pointLevelPtr_) && !pointLevelPtr_)
    {
        IOobject rio(io, "pointLevel");
        rio.readOpt(IOobject::NO_READ);
        pointLevelPtr_.reset
        (
            new labelIOList(rio, labelList(mesh.nPoints(), Zero))
        );
    }
    // Take over IOobject settings
    if (pointLevelPtr_)
    {
        pointLevelPtr_->IOobjectOption::operator=(io);
    }

    if (returnReduceOr(level0EdgePtr_))
    {
        // Get master length
        scalar masterLen = (UPstream::master() ? level0EdgePtr_().value() : 0);
        Pstream::broadcast(masterLen);
        if (!level0EdgePtr_)
        {
            IOobject rio(io, "level0Edge");
            rio.readOpt(IOobject::NO_READ);
            level0EdgePtr_.reset
            (
                new uniformDimensionedScalarField(rio, dimLength, masterLen)
            );
        }
    }
    // Take over IOobject settings
    if (level0EdgePtr_)
    {
        level0EdgePtr_->IOobjectOption::operator=(io);
    }

    if (returnReduceOr(refHistoryPtr_) && !refHistoryPtr_)
    {
        IOobject rio(io, "refinementHistory");
        rio.readOpt(IOobject::NO_READ);
        refHistoryPtr_.reset(new refinementHistory(rio, mesh.nCells(), true));
    }
    // Take over IOobject settings
    if (refHistoryPtr_)
    {
        refHistoryPtr_->IOobjectOption::operator=(io);
    }
}


void Foam::hexRef8Data::updateMesh(const mapPolyMesh& map)
{
    // Sanity check
    if
    (
         (cellLevelPtr_ && cellLevelPtr_().size() != map.nOldCells())
      || (pointLevelPtr_ && pointLevelPtr_().size() != map.nOldPoints())
    )
    {
        cellLevelPtr_.clear();
        pointLevelPtr_.clear();
        level0EdgePtr_.clear();
        refHistoryPtr_.clear();
        return;
    }


    if (cellLevelPtr_)
    {
        const labelList& cellMap = map.cellMap();
        labelList& cellLevel = cellLevelPtr_();

        labelList newCellLevel(cellMap.size());
        forAll(cellMap, newCelli)
        {
            label oldCelli = cellMap[newCelli];

            if (oldCelli == -1)
            {
                newCellLevel[newCelli] = 0;
            }
            else
            {
                newCellLevel[newCelli] = cellLevel[oldCelli];
            }
        }
        cellLevel.transfer(newCellLevel);
        cellLevelPtr_().instance() = map.mesh().facesInstance();
    }
    if (pointLevelPtr_)
    {
        const labelList& pointMap = map.pointMap();
        labelList& pointLevel = pointLevelPtr_();

        labelList newPointLevel(pointMap.size());
        forAll(pointMap, newPointi)
        {
            label oldPointi = pointMap[newPointi];

            if (oldPointi == -1)
            {
                newPointLevel[newPointi] = 0;
            }
            else
            {
                newPointLevel[newPointi] = pointLevel[oldPointi];
            }
        }
        pointLevel.transfer(newPointLevel);
        pointLevelPtr_().instance() = map.mesh().facesInstance();
    }


    if (refHistoryPtr_ && refHistoryPtr_().active())
    {
        refHistoryPtr_().updateMesh(map);
        refHistoryPtr_().instance() = map.mesh().facesInstance();
    }
}


void Foam::hexRef8Data::distribute(const mapDistributePolyMesh& map)
{
    if (cellLevelPtr_)
    {
        map.cellMap().distribute(*cellLevelPtr_);
    }
    if (pointLevelPtr_)
    {
        map.pointMap().distribute(*pointLevelPtr_);
    }

    // No need to distribute the level0Edge

    if (refHistoryPtr_ && refHistoryPtr_().active())
    {
        refHistoryPtr_().distribute(map);
    }
}


bool Foam::hexRef8Data::write() const
{
    bool ok = true;
    if (cellLevelPtr_)
    {
        ok = ok && cellLevelPtr_().write();
    }
    if (pointLevelPtr_)
    {
        ok = ok && pointLevelPtr_().write();
    }
    if (level0EdgePtr_)
    {
        ok = ok && level0EdgePtr_().write();
    }
    if (refHistoryPtr_)
    {
        ok = ok && refHistoryPtr_().write();
    }
    return ok;
}


// ************************************************************************* //

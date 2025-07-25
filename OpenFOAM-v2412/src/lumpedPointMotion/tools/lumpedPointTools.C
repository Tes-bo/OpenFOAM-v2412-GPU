/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "lumpedPointTools.H"
#include "IFstream.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "points0MotionSolver.H"
#include "lumpedPointDisplacementPointPatchVectorField.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class GeoFieldType>
static autoPtr<GeoFieldType> loadPointField
(
    const pointMesh::Mesh& mesh,
    const IOobject* io
)
{
    if (io && io->isHeaderClass<GeoFieldType>())
    {
        Info<< "Reading " << io->headerClassName()
            << ' ' << io->name() << endl;

        return autoPtr<GeoFieldType>::New
        (
            IOobject
            (
                io->name(),
                io->instance(),
                io->local(),
                io->db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                io->registerObject()
            ),
            mesh
        );
    }

    return nullptr;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::lumpedPointStateTuple>
Foam::lumpedPointTools::lumpedPointStates
(
    const dictionary& dict,
    quaternion::eulerOrder rotOrder,
    bool degrees
)
{
    quaternion::eulerOrderNames.readIfPresent
    (
        "rotationOrder",
        dict,
        rotOrder
    );

    dict.readIfPresent("degrees", degrees);

    Info<<"Reading states\n";
    List<dictionary> entries(dict.lookup("response"));

    DynamicList<Tuple2<scalar, lumpedPointState>> states(entries.size());

    for (const dictionary& subDict : entries)
    {
        states.append
        (
            lumpedPointStateTuple
            (
                subDict.get<scalar>("time"),
                lumpedPointState(subDict)
            )
        );
    }

    return states.shrink();
}


Foam::List<Foam::lumpedPointStateTuple>
Foam::lumpedPointTools::lumpedPointStates
(
    Istream& is,
    quaternion::eulerOrder rotOrder,
    bool degrees
)
{
    dictionary dict(is);
    return lumpedPointStates(dict, rotOrder, degrees);
}


Foam::List<Foam::lumpedPointStateTuple>
Foam::lumpedPointTools::lumpedPointStates
(
    const fileName& file,
    quaternion::eulerOrder rotOrder,
    bool degrees
)
{
    IFstream is(file);
    return lumpedPointStates(is, rotOrder, degrees);
}


Foam::pointIOField
Foam::lumpedPointTools::points0Field(const polyMesh& mesh)
{
    return pointIOField(points0MotionSolver::points0IO(mesh));
}


Foam::labelList
Foam::lumpedPointTools::lumpedPointPatchList(const pointVectorField& pvf)
{
    return lumpedPointDisplacementPointPatchVectorField::patchIds(pvf);
}


Foam::labelList Foam::lumpedPointTools::lumpedPointPatchList
(
    const polyMesh& mesh
)
{
    IOobjectList objects0(mesh, "0");

    pointMesh pMesh(mesh);

    autoPtr<pointVectorField> displacePtr =
        loadPointField<pointVectorField>
        (
            pMesh,
            objects0.findObject("pointDisplacement")
        );

    if (!displacePtr)
    {
        Info<< "No valid pointDisplacement" << endl;
        return labelList();
    }

    return lumpedPointPatchList(*displacePtr);
}


Foam::label
Foam::lumpedPointTools::setPatchControls
(
    const pointVectorField& pvf,
    const pointField& points0
)
{
    return
        lumpedPointDisplacementPointPatchVectorField::setPatchControls
        (
            pvf,
            points0
        );
}


Foam::label Foam::lumpedPointTools::setPatchControls
(
    const fvMesh& mesh,
    const pointField& points0
)
{
    IOobjectList objects0(mesh, "0");

    pointMesh pMesh(mesh);

    autoPtr<pointVectorField> displacePtr =
        loadPointField<pointVectorField>
        (
            pMesh,
            objects0.findObject("pointDisplacement")
        );

    if (!displacePtr)
    {
        Info<< "No valid pointDisplacement" << endl;
        return 0;
    }

    return setPatchControls(*displacePtr, points0);
}


Foam::label Foam::lumpedPointTools::setPatchControls
(
    const fvMesh& mesh
)
{
    pointIOField points0(points0Field(mesh));

    return setPatchControls(mesh, points0);
}


Foam::label Foam::lumpedPointTools::setInterpolators
(
    const pointVectorField& pvf,
    const pointField& points0
)
{
    return
        lumpedPointDisplacementPointPatchVectorField::setInterpolators
        (
            pvf,
            points0
        );
}


Foam::label Foam::lumpedPointTools::setInterpolators
(
    const fvMesh& mesh,
    const pointField& points0
)
{
    IOobjectList objects0(mesh, "0");

    pointMesh pMesh(mesh);

    autoPtr<pointVectorField> displacePtr =
        loadPointField<pointVectorField>
        (
            pMesh,
            objects0.findObject("pointDisplacement")
        );

    if (!displacePtr)
    {
        Info<< "No valid pointDisplacement" << endl;
        return 0;
    }

    return setInterpolators(*displacePtr, points0);
}


Foam::label Foam::lumpedPointTools::setInterpolators
(
    const fvMesh& mesh
)
{
    pointIOField points0(points0Field(mesh));

    return setInterpolators(mesh, points0);
}


// ************************************************************************* //

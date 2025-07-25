/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

Notes
    Renumber using Zoltan

    Zoltan install:
    - in your ~/.bashrc:
            export ZOLTAN_ARCH_DIR=\
                $WM_THIRD_PARTY_DIR/platforms/linux64Gcc/Zoltan_XXX
    - unpack into $WM_THIRD_PARTY_DIR
    - cd Zoltan_XXX
    - mkdir build
    - cd build
    - export CCFLAGS="-fPIC"
    - export CXXFLAGS="-fPIC"
    - export CFLAGS="-fPIC"
    - export LDFLAGS="-shared"
    - ../configure \
        --prefix=$ZOLTAN_ARCH_DIR \
        --with-ccflags=-fPIC --with-cxxflags=-fPIC --with-ldflags=-shared

\*---------------------------------------------------------------------------*/

#include "zoltanRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "labelIOList.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "CStringList.H"
#include "uint.H"
#include <algorithm>
#include <numeric>

// Include MPI without any C++ bindings
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif

#pragma GCC diagnostic ignored "-Wold-style-cast"
#include "zoltan.h"
#include <mpi.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zoltanRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        zoltanRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// [ZOLTAN_NUM_OBJ_FN]
// The number of graph vertices (locally)
static int get_number_of_vertices(void *data, int *ierr)
{
    *ierr = ZOLTAN_OK;
    const auto& mesh = *static_cast<const Foam::polyMesh*>(data);

    return mesh.nCells();
}


// [ZOLTAN_OBJ_LIST_FN]
// The ids for the graph vertices
static void get_vertex_list
(
    void *data,
    int sizeGID,                // (currently unused)
    int sizeLID,                // (currently unused)
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int wgt_dim,                // (currently unused)
    float *obj_wgts,            // (currently unused)
    int *ierr
)
{
    *ierr = ZOLTAN_OK;
    const auto& mesh = *static_cast<const Foam::polyMesh*>(data);

    const auto nCells = mesh.nCells();

    // Offset for globally unique IDs
    const auto myProci = Foam::UPstream::myProcNo();
    const auto myProcOffset =
        mesh.globalData().globalMeshCellAddr().localStart(myProci);

    // Global indices
    std::iota(global_ids, (global_ids + nCells), ZOLTAN_ID_TYPE(myProcOffset));

    // Local indices
    std::iota(local_ids, (local_ids + nCells), ZOLTAN_ID_TYPE(0));

    // No weights?
    // Zoltan will assume equally weighted vertices.

    // if (obj_wgts && wgt_dim > 0)
    // {
    //     std::fill_n(obj_wgts, wgt_dim, float(1));
    // }
}


// [ZOLTAN_NUM_EDGES_MULTI_FN]
// query function returns the number of edges in the communication
// graph of the application for each object in a list of objects.
// That is, for each object in the global_ids/local_ids arrays, the
// number of objects with which the given object must share
// information is returned.

static void get_num_edges_list
(
    void *data,
    int sizeGID,
    int sizeLID,
    int num_obj,                // Number of graph vertices (mesh cells)
    ZOLTAN_ID_PTR global_ids,   // (currently unused)
    ZOLTAN_ID_PTR local_ids,    // rank-local vertex id (mesh cell id)
    int *numEdges,
    int *ierr
)
{
    *ierr = ZOLTAN_OK;
    const auto& mesh = *static_cast<const Foam::polyMesh*>(data);

    if ((sizeGID != 1) || (sizeLID != 1) || (num_obj != mesh.nCells()))
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    const Foam::label nCells = num_obj;

    for (Foam::label i=0; i < nCells; ++i)
    {
        const Foam::label celli = local_ids[i];
        int numNbr = 0;

        for (const auto facei : mesh.cells()[celli])
        {
            if (mesh.isInternalFace(facei))
            {
                ++numNbr;
            }
            // TBD: check coupled etc
        }

        numEdges[i] = numNbr;
    }
}


// [ZOLTAN_EDGE_LIST_MULTI_FN]
static void get_edge_list
(
    void *data,
    int sizeGID,
    int sizeLID,
    int num_obj,
    ZOLTAN_ID_PTR global_ids,   // (currently unused)
    ZOLTAN_ID_PTR local_ids,    // rank-local vertex id (mesh cell id)
    int *num_edges,
    ZOLTAN_ID_PTR nborGID,
    int *nborProc,
    int wgt_dim,
    float *ewgts,
    int *ierr
)
{
    *ierr = ZOLTAN_OK;
    const auto& mesh = *static_cast<const Foam::polyMesh*>(data);

    if
    (
        (sizeGID != 1)
     || (sizeLID != 1)
     || (num_obj != mesh.nCells())
     || (wgt_dim != 1)
    )
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    // Offset for globally unique IDs
    const auto myProci = Foam::UPstream::myProcNo();
    const auto myProcOffset =
        mesh.globalData().globalMeshCellAddr().localStart(myProci);

    auto* nextNbor = nborGID;
    auto* nextProc = nborProc;
    auto* nextWgt = ewgts;

    const Foam::label nCells = num_obj;

    for (Foam::label i=0; i < nCells; ++i)
    {
        const Foam::label celli = local_ids[i];
        int numNbr = 0;

        for (const auto facei : mesh.cells()[celli])
        {
            if (mesh.isInternalFace(facei))
            {
                Foam::label nbr = mesh.faceOwner()[facei];
                if (nbr == celli)
                {
                    nbr = mesh.faceNeighbour()[facei];
                }

                *nextNbor++ = (nbr + myProcOffset); // global id
                *nextProc++ = myProci;              // rank-local connection
                *nextWgt++ = 1.0;

                ++numNbr;
            }
        }

        // Sanity check
        if (numNbr != num_edges[i])
        {
            *ierr = ZOLTAN_FATAL;
            return;
        }
    }
}


// [ZOLTAN_NUM_GEOM_FN]
// The dimensionality of the mesh geometry (2D,3D, etc)
static int get_mesh_dim(void *data, int *ierr)
{
    *ierr = ZOLTAN_OK;
    const auto& mesh = *static_cast<const Foam::polyMesh*>(data);

    return mesh.nSolutionD();
}


// [ZOLTAN_GEOM_MULTI_FN]
// The geometric location of the graph vertices (mesh cellCentres)
static void get_geom_list
(
    void *data,
    int num_gid_entries,
    int num_lid_entries,
    int num_obj,
    ZOLTAN_ID_PTR global_ids,
    ZOLTAN_ID_PTR local_ids,
    int num_dim,                // dimensionality (2|3)
    double *geom_vec,           // [out] cellCentres
    int *ierr
)
{
    *ierr = ZOLTAN_OK;
    const Foam::polyMesh& mesh = *static_cast<const Foam::polyMesh*>(data);

    if
    (
        (num_gid_entries != 1)
     || (num_lid_entries != 1)
     || (num_obj != mesh.nCells())
     || (num_dim != mesh.nSolutionD())
    )
    {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    const Foam::pointField& cc = mesh.cellCentres();

    if (num_dim == 3)
    {
        // Fast path for 3D - assumes that local_ids are still ordered
        std::copy_n
        (
            reinterpret_cast<const Foam::point::cmptType*>(cc.cdata()),
            (3*num_obj),
            geom_vec
        );
    }
    else
    {
        // [out] cellCentres
        double* p = geom_vec;

        const Foam::Vector<Foam::label>& sol = mesh.solutionD();

        const Foam::label nCells = num_obj;

        // Assumes that local_ids are still ordered
        for (Foam::label celli = 0; celli < nCells; ++celli)
        {
            const Foam::point& pt = cc[celli];

            for
            (
                Foam::direction cmpt = 0;
                cmpt < Foam::vector::nComponents;
                ++cmpt
            )
            {
                if (sol[cmpt] == 1)
                {
                    *p++ = pt[cmpt];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoltanRenumber::zoltanRenumber(const dictionary& dict)
:
    renumberMethod(dict),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::zoltanRenumber::renumber
(
    const polyMesh& pMesh
) const
{
    // Zoltan_Initialize will trigger MPI_Init() if not already done.
    // - use UPstream::initNull() so that OpenFOAM also knows about MPI

    UPstream::initNull();

    CStringList args({"zoltan-renumber"});

    float ver;
    int rc = Zoltan_Initialize(args.size(), args.strings(), &ver);

    if (rc != ZOLTAN_OK)
    {
        FatalErrorInFunction
            << "Failed initialising Zoltan" << exit(FatalError);
    }

    struct Zoltan_Struct *zz = Zoltan_Create(MPI_COMM_WORLD);

    // const bool verbose = coeffsDict_.getOrDefault(verbose, false);
    const bool verbose = true;

    // Default order method
    Zoltan_Set_Param(zz, "ORDER_METHOD", "LOCAL_HSFC");

    if (false)
    {
        Info<< typeName << " : default ORDER_METHOD = LOCAL_HSFC" << nl
            << typeName << " : default ORDER_TYPE = LOCAL" << nl;
    }

    for (const entry& dEntry : coeffsDict_)
    {
        const word& key = dEntry.keyword();

        // Internal keywords
        if
        (
            key == "method"
         || key == "verbose"
        )
        {
            continue;
        }

        if (!dEntry.isDict())
        {
            const word value(dEntry.get<word>());

            if (verbose)
            {
                Info<< typeName
                    << " : setting parameter "
                    << key << " = " << value << nl;
            }

            Zoltan_Set_Param(zz, key.c_str(), value.c_str());
        }
    }

    // Always use rank LOCAL ordering
    Zoltan_Set_Param(zz, "ORDER_TYPE", "LOCAL");


    // Set callbacks
    polyMesh& mesh = const_cast<polyMesh&>(pMesh);

    // Callbacks for graph vertex IDs
    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &mesh);
    Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &mesh);

    // Callbacks for geometry
    Zoltan_Set_Num_Geom_Fn(zz, get_mesh_dim, &mesh);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geom_list, &mesh);

    // Callbacks for connectivity
    Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges_list, &mesh);
    Zoltan_Set_Edge_List_Multi_Fn(zz, get_edge_list, &mesh);


    const auto nCells = mesh.nCells();

    List<ZOLTAN_ID_TYPE> globalIds(nCells);
    List<ZOLTAN_ID_TYPE> oldToNew(nCells);

    // Offset for globally unique IDs
    const label myProci = UPstream::myProcNo();
    const label myProcOffset =
        mesh.globalData().globalMeshCellAddr().localStart(myProci);


    // Global indices
    std::iota(globalIds.begin(), globalIds.end(), ZOLTAN_ID_TYPE(myProcOffset));

    int err = Zoltan_Order
    (
        zz,
        1,                          //int num_gid_entries,
        nCells,                     //int num_obj,
        globalIds.begin(),
        oldToNew.begin()
    );

    if (err != ZOLTAN_OK)
    {
        FatalErrorInFunction
            << "Failed Zoltan_Order" << exit(FatalError);
    }

    Zoltan_Destroy(&zz);


    // From global number to rank-local numbering
    globalIds.clear();
    labelList order(oldToNew.size());

    // From global to local (without checks)
    std::transform
    (
        oldToNew.begin(),
        oldToNew.end(),
        order.begin(),
        [=](const ZOLTAN_ID_TYPE id) -> label { return (id - myProcOffset); }
    );

    return order;
}


// ************************************************************************* //

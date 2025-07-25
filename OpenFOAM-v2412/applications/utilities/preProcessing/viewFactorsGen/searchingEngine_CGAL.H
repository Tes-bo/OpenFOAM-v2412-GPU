Random rndGen(653213);

// Determine mesh bounding boxes:
List<treeBoundBox> meshBb
(
    1,
    treeBoundBox(coarseMesh.points()).extend(rndGen, 1e-3)
);

// Dummy bounds dictionary
dictionary dict;
dict.add("bounds", meshBb);
dict.add
(
    "distributionType",
    distributedTriSurfaceMesh::distributionTypeNames_
    [
        distributedTriSurfaceMesh::FROZEN
    ]
);
dict.add("mergeDistance", SMALL);

labelList triSurfaceToAgglom(5*nFineFaces);

const triSurface localSurface = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches
);


// CGAL surface

distributedTriSurfaceMesh surfacesMesh
(
    IOobject
    (
        "wallSurface.stl",
        runTime.constant(),     // directory
        "triSurface",           // instance
        runTime,                // registry
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    localSurface,
    dict
);

triSurfaceToAgglom.resize(surfacesMesh.size());

surfacesMesh.setField(triSurfaceToAgglom);

std::vector<Triangle> triangles;

{
    const auto& surf = static_cast<const triSurface&>(surfacesMesh);

    triangles.reserve(surf.size());

    const pointField& pts = surf.points();

    for (const auto& f : surf)
    {
        const point& a = pts[f.a()];
        const point& b = pts[f.b()];
        const point& c = pts[f.c()];

        triangles.emplace_back
        (
            Point(a.x(), a.y(), a.z()),
            Point(b.x(), b.y(), b.z()),
            Point(c.x(), c.y(), c.z())
        );

        if (triangles.back().is_degenerate())
        {
            std::cout << triangles.back() << std::endl;
        }
    }
}


// constructs AABB tree
Tree tree(triangles.begin(), triangles.end());

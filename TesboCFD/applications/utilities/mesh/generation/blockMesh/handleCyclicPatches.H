// Detect any cyclic patches and force re-ordering of the faces
{
    bool hasCyclic = false;
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (isA<cyclicPolyPatch>(pp))
        {
            hasCyclic = true;
            break;
        }
    }

    if (hasCyclic)
    {
        Info<< nl << "Detected cyclic patches; ordering boundary faces" << endl;

        const word oldInstance = mesh.instance();
        polyTopoChange meshMod(mesh);
        meshMod.changeMesh(mesh, false);
        mesh.setInstance(oldInstance);
    }
}

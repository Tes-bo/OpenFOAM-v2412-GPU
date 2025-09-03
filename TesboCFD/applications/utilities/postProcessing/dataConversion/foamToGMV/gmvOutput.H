for (const word& fieldType : fieldTypes)
{
    const wordList fieldNames = objects.sortedNames(fieldType);

    if (fieldType == "volScalarField")
    {
        gmvFile << "variable" << nl;
    }

    for (const word& fieldName : fieldNames)
    {
        IOobject fieldObject
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldType == "volScalarField")
        {
            volScalarField fld(fieldObject, mesh);
            gmvFile << fieldName << " 0" << nl;
            for (label indx=0; indx<mesh.nCells(); ++indx)
            {
                gmvFile << fld[indx] << " ";
            }
            gmvFile << nl;
        }

        if (fieldType == "volVectorField")
        {
            if (fieldName == vComp)
            {
                volVectorField fld(fieldObject, mesh);
                gmvFile << "velocity 0" << nl;
                for (label indx=0;indx<mesh.nCells();indx++)
                {
                    gmvFile << fld[indx].x() << " ";
                }
                for (label indx=0;indx<mesh.nCells();indx++)
                {
                    gmvFile << fld[indx].y() << " ";
                }
                for (label indx=0;indx<mesh.nCells();indx++)
                {
                    gmvFile << fld[indx].z() << " ";
                }
                gmvFile << nl;
            }
        }

        if (fieldType == "surfaceScalarField")
        {
            // ...
        }

    }

    if (fieldType == cloud::prefix)
    {
        IOobject positionsHeader
        (
            "positions",
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::NO_READ
        );

        IOobject coordinatesHeader
        (
            "coordinates",
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::NO_READ
        );

        if
        (
            positionsHeader.typeHeaderOk<IOPosition<Cloud<passiveParticle>>>
            (
                false
            )
         || coordinatesHeader.typeHeaderOk<IOPosition<Cloud<passiveParticle>>>
            (
                false
            )
        )
        {
            Cloud<passiveParticle> particles(mesh, cloud::defaultName);

            IOobjectList objects(mesh, runTime.timeName(), cloud::prefix);

            wordList lagrangianScalarNames(objects.sortedNames("scalarField"));
            wordList lagrangianVectorNames(objects.sortedNames("vectorField"));

            if (particles.size())
            {
                #include "gmvOutputLagrangian.H"
            }
        }
    }

    if (fieldType == "volScalarField")
    {
        gmvFile << "endvars" << nl;
    }
}

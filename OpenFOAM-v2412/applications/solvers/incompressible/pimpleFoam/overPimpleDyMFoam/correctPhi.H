if (mesh.changing())
{
    volVectorField::Boundary& bfld = U.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].initEvaluate();
        }
    }

    surfaceScalarField::Boundary& phiBfld = phi.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].evaluate();

            phiBfld[patchi] = bfld[patchi] & mesh.Sf().boundaryField()[patchi];
        }
    }
}
    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        fvPatchFieldBase::zeroGradientType()
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(p.dimensions(), Zero),
        pcorrTypes
    );

{
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    const cellCellStencilObject& overlap = Stencil::New(mesh);
    const labelList& cellTypes = overlap.cellTypes();
    const labelIOList& zoneIDs = overlap.zoneID();

    while (pimple.correctNonOrthogonal())
    {
        label nZones = gMax(zoneIDs)+1;

        //label refCellI2 = -1;
        labelList refCells(nZones, -1);
        labelList refZones(nZones, -1);

        forAll(zoneIDs, cellI)
        {
            label zoneId = zoneIDs[cellI];
            if
            (
                refCells[zoneId] == -1
             && cellTypes[cellI] == cellCellStencil::CALCULATED
             && refZones[zoneId] == -1
            )
            {
                refCells[zoneId] = cellI;
                refZones[zoneId] = zoneId;
            }
        }

        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAUf, pcorr) == fvc::div(phi)
        );


        // Only set reference for cells that are CALCULATED
        {
            DynamicList<label> validCells(refCells.size());
            forAll(refCells, zoneId)
            {
                if (refCells[zoneId] != -1)
                {
                    validCells.push_back(refCells[zoneId]);
                }
            }

            pcorrEqn.setReferences
            (
                validCells,
                scalar(0),
                true
            );
        }

        const dictionary& d = mesh.solver
        (
            pcorr.select
            (
                pimple.finalInnerIter()
            )
        );
        //Bypass virtual layer
        mesh.fvMesh::solve(pcorrEqn, d);

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }
}

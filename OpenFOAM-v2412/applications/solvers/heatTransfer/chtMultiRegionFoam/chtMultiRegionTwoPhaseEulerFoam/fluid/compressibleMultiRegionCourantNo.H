    scalar CoNum = -GREAT;

    forAll(fluidRegions, regioni)
    {
        if (fluidRegions[regioni].nInternalFaces())
        {
            const surfaceScalarField& phi =
                phaseSystemFluid[regioni].phi();

            scalarField sumPhi
            (
                fvc::surfaceSum(mag(phi))().primitiveField()
            );

            const surfaceScalarField& phi1 =
                phaseSystemFluid[regioni].phase1().phiRef();

            const surfaceScalarField& phi2 =
                phaseSystemFluid[regioni].phase2().phiRef();

            sumPhi = max
            (
                sumPhi,
                fvc::surfaceSum(mag(phi1))().primitiveField()
            );

            sumPhi = max
            (
                sumPhi,
                fvc::surfaceSum(mag(phi2))().primitiveField()
            );


            CoNum =
                0.5*gMax
                (
                    sumPhi/fluidRegions[regioni].V().field()
                )*runTime.deltaTValue();

            scalar UrCoNum = 0.5*gMax
            (
                fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()
              / fluidRegions[regioni].V().field()
            )*runTime.deltaTValue(),

            CoNum = max(UrCoNum, CoNum);
        }
    }

    Info<< "Courant Number max: " << CoNum << endl;

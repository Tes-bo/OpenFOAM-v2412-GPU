    dimensionedScalar StCorr("StCorr", dimless, 1.0);

    if (ign.igniting())
    {
        // Calculate volume of ignition kernel
        dimensionedScalar Vk("Vk", dimVolume, gSum(c*mesh.V().field()));
        dimensionedScalar Ak("Ak", dimArea, Zero);

        if (Vk.value() > SMALL)
        {
            // Calculate kernel area from its volume
            // and the dimensionality of the case

            switch(mesh.nGeometricD())
            {
                case 3:
                {
                    // Assume it is part-spherical
                    scalar sphereFraction
                    (
                        combustionProperties.get<scalar>
                        (
                            "ignitionSphereFraction"
                        )
                    );

                    Ak = sphereFraction*4.0*constant::mathematical::pi
                       *pow
                        (
                            3.0*Vk
                           /(sphereFraction*4.0*constant::mathematical::pi),
                            2.0/3.0
                        );
                }
                break;

                case 2:
                {
                    // Assume it is part-circular
                    dimensionedScalar thickness
                    (
                        "ignitionThickness", dimLength, combustionProperties
                    );

                    scalar circleFraction
                    (
                        combustionProperties.get<scalar>
                        (
                            "ignitionCircleFraction"
                        )
                    );

                    Ak = circleFraction*constant::mathematical::pi*thickness
                       *sqrt
                        (
                            4.0*Vk
                           /(
                               circleFraction
                              *thickness
                              *constant::mathematical::pi
                            )
                        );
                }
                break;

                case 1:
                    // Assume it is plane or two planes
                    Ak = dimensionedScalar
                    (
                        "ignitionKernelArea", dimArea, combustionProperties
                    );
                break;
            }

            // Calculate kernel area from b field consistent with the
            // discretisation of the b equation.
            const volScalarField mgb
            (
                fvc::div(nf, b, "div(phiSt,b)") - b*fvc::div(nf) + dMgb
            );
            dimensionedScalar AkEst = gSum(mgb*mesh.V().field());

            StCorr.value() = clamp((Ak/AkEst).value(), 1, 10);

            Info<< "StCorr = " << StCorr.value() << endl;
        }
    }

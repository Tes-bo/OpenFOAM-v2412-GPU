// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Process some "standard" list options

#ifndef setRootCaseListOutput_H
#define setRootCaseListOutput_H

{
    bool listOptions = false;

    if (args.found("listSwitches"))
    {
        debug::listSwitches(args.found("listUnsetSwitches"));
        listOptions = true;
    }

    if (args.found("listRegisteredSwitches"))
    {
        debug::listRegisteredSwitches(args.found("listUnsetSwitches"));
        listOptions = true;
    }

    #ifdef Foam_fvPatchField_H
    if (args.found("listScalarBCs"))
    {
        Info<< "scalarBCs"
            << fvPatchField<Foam::scalar>::
               dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        listOptions = true;
    }

    if (args.found("listVectorBCs"))
    {
        Info<< "vectorBCs"
            <<  fvPatchField<Foam::vector>::
                dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    #ifdef Foam_functionObject_H
    if (args.found("listFunctionObjects"))
    {
        Info<< "functionObjects"
            << functionObject::dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    #ifdef Foam_fvOption_H
    if (args.found("listFvOptions"))
    {
        Info<< "fvOptions"
            << fv::option::dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    #if defined(Foam_turbulentTransportModel_H)
    if (args.found("listTurbulenceModels"))
    {
        Info<< "Turbulence models" << nl
            << incompressible::turbulenceModel::
               dictionaryConstructorTablePtr_->sortedToc()
            << nl << endl;

        Info<< "RAS models"
            << incompressible::RASModel::
               dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "LES models"
            << incompressible::LESModel::
               dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        listOptions = true;
    }
    #elif defined(Foam_turbulentFluidThermoModel_H)
    if (args.found("listTurbulenceModels"))
    {
        Info<< "Turbulence models" << nl
            << compressible::turbulenceModel::
               dictionaryConstructorTablePtr_->sortedToc()
            << nl << endl;

        Info<< "RAS models"
            << compressible::RASModel::
               dictionaryConstructorTablePtr_->sortedToc()
            << endl;

        Info<< "LES models"
            << compressible::LESModel::
               dictionaryConstructorTablePtr_->sortedToc()
            << endl;
        listOptions = true;
    }
    #endif

    if (listOptions)
    {
        exit(0);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

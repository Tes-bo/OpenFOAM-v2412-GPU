    wordList pcorrTypes
    (
        p.boundaryField().size(),
        fvPatchFieldBase::zeroGradientType()
    );

    for (label i=0; i<p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

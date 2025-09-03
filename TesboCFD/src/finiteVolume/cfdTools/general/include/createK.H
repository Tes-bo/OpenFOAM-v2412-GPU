Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

if (U.nOldTimes())
{
    volVectorField* Uold = &U.oldTime();
    volScalarField* Kold = &K.oldTime();
    *Kold == 0.5*magSqr(*Uold);

    while (Uold->nOldTimes())
    {
        Uold = &Uold->oldTime();
        Kold = &Kold->oldTime();
        *Kold == 0.5*magSqr(*Uold);
    }
}

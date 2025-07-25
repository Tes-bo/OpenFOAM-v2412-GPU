Info<< "Reading field Cs" << endl;
areaScalarField Cs
(
    IOobject
    (
        "Cs",
        runTime.timeName(),
        aMesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);

dimensionedScalar Cs0("Cs0", dimMass/dimArea, 1.0);

const areaVectorField& R = aMesh.areaCentres();

Cs = Cs0*(1.0 + R.component(vector::X)/mag(R));

dimensionedScalar Ds("Ds", dimViscosity, 1.0);

areaVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        aMesh.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    aMesh,
    dimensionedVector(dimVelocity, Zero)
);

dimensioned<scalar> Uinf("Uinf", dimVelocity, 1.0);

forAll(Us, faceI)
{
    Us[faceI].x() =
        Uinf.value()*(0.25*(3.0 + sqr(R[faceI].x()/mag(R[faceI]))) - 1.0);

    Us[faceI].y() =
        Uinf.value()*0.25*R[faceI].x()*R[faceI].y()/sqr(mag(R[faceI]));

    Us[faceI].z() =
        Uinf.value()*0.25*R[faceI].x()*R[faceI].z()/sqr(mag(R[faceI]));
}
Us.boundaryFieldRef().evaluateCoupled<coupledFaPatch>();

Us -= aMesh.faceAreaNormals()*(aMesh.faceAreaNormals() & Us);


edgeScalarField phis
(
    IOobject
    (
        "phis",
        runTime.timeName(),
        aMesh.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearEdgeInterpolate(Us) & aMesh.Le()
);

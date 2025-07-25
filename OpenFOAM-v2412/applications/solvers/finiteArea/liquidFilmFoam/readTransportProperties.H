IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar mug("mug", dimViscosity, transportProperties);
dimensionedScalar mul("mul", dimViscosity, transportProperties);

dimensionedScalar rhog("rhog", dimDensity, transportProperties);
dimensionedScalar rhol("rhol", dimDensity, transportProperties);

dimensionedScalar sigma("sigma", dimMass/sqr(dimTime), transportProperties);

dimensionedScalar h0("h0", dimLength, transportProperties);

{
    // Stabilisation of friction factor calculation
    // Friction factor is defined with standard gravity
    frictionFactor.primitiveFieldRef() =
        mag(2*9.81*sqr(manningField.primitiveField())/
        cbrt(mag(h.primitiveField()) + 1e-7));
}

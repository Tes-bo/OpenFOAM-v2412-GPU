{
    scalar CoNumSigma = max
    (
        sqrt
        (
            2*M_PI*sigma*sqr(aMesh.edgeInterpolation::deltaCoeffs())
            *mag(aMesh.edgeInterpolation::deltaCoeffs())
            /rhol
        )
    ).value()*runTime.deltaTValue();

    Info<< "Max Capillary Courant Number = " << CoNumSigma << '\n' << endl;
}

## Implicit flux conservation test

The implicit flux conservation is switched on in the field (here T) overset BC.

The keyword is used as:

    free
    {
        type            overset;
        fluxCorrection  true;
        value           uniform 300;
    }

A special net fringe flux is output in this tutorial.

Switch `fluxCorrection` to `false` to evaluate the flux correction.


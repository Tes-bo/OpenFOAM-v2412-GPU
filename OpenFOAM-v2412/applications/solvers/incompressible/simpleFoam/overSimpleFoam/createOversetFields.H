//- Overset specific

// Add solver-specific interpolations
{
    wordHashSet& nonInt =
        const_cast<wordHashSet&>(Stencil::New(mesh).nonInterpolatedFields());

    nonInt.insert("HbyA");
    nonInt.insert("grad(p)");
    nonInt.insert("surfaceIntegrate(phi)");
    nonInt.insert("surfaceIntegrate(phiHbyA)");
    nonInt.insert("cellMask");
    nonInt.insert("cellDisplacement");
    nonInt.insert("interpolatedCells");
    nonInt.insert("cellInterpolationWeight");
}

// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"

#include "createInterpolatedCells.H"

bool adjustFringe
(
    simple.dict().getOrDefault("oversetAdjustPhi", false)
);

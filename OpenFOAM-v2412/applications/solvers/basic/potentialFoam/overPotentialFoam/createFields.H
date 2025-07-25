#include "../createFields.H"

// Add solver-specific interpolations
{
    wordHashSet& nonInt =
        const_cast<wordHashSet&>(Stencil::New(mesh).nonInterpolatedFields());

    nonInt.insert("cellMask");
    nonInt.insert("interpolatedCells");

}

// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"

// Create bool field with interpolated cells
#include "createInterpolatedCells.H"

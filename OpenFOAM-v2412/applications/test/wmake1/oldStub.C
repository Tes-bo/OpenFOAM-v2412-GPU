// Some test code

#include "IOstreams.H"

namespace Foam
{
    void printTest()
    {
        Info<< nl << "Using old stub" << nl;

        #if OPENFOAM
        Info<< "OPENFOAM=" << OPENFOAM << nl;
        #else
        Info<< "OPENFOAM is undefined" << nl;
        #endif
    }
}

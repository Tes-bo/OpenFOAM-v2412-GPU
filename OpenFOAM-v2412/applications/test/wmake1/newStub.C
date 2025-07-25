// Some test code

#include "foamVersion.H"
#include "IOstreams.H"

namespace Foam
{
    void printTest()
    {
        Info<< nl;
        foamVersion::printBuildInfo(Info.stdStream());
    }
}

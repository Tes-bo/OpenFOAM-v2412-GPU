#include "argList.H"
#include "point.H"
#include "triangle.H"
#include "tetrahedron.H"
#include "IOstreams.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    triangle<point, point> tri
    (
        vector(0, 0, 0),
        vector(1, 0, 0),
        vector(1, 1, 0)
    );

    Info<< "triangle: " << tri << nl
        << "  vecA: " << tri.vecA() << nl
        << "  vecB: " << tri.vecB() << nl
        << "  vecC: " << tri.vecC() << nl
        << endl;

    Info<< "tri circumCentre = " << tri.circumCentre() << endl;
    Info<< "tri circumRadius = " << tri.circumRadius() << endl;

    tetrahedron<point, point> tet
    (
        vector(1, 0, 0),
        vector(0, 1, 0),
        vector(0, 0, 1),
        vector(0.5773502, 0.5773502, 0.5773502)
    );

    Info<< "tet circumCentre = " << tet.circumCentre() << endl;
    Info<< "tet circumRadius = " << tet.circumRadius() << endl;

    InfoErr<< "Enter four points: " << endl;

    vector a(Sin);
    vector b(Sin);
    vector c(Sin);
    vector d(Sin);

    Info<< "tet circumRadius = "
        << tetrahedron<point, point>(a, b, c, d).circumRadius() << endl;

    Info<< "\nEnd\n";
    return 0;
}

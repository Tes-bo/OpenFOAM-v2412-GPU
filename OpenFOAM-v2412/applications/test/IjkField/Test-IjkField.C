/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Functionality of IjkField

\*---------------------------------------------------------------------------*/

#include "point.H"
#include "IjkField.H"
#include "IOstreams.H"

using namespace Foam;

template<class T>
Ostream& print(const IjkField<T>& fld)
{
    Info<< static_cast<const Field<T>&>(fld).size()
        << " addr:" << name(fld.cdata()) << ' ' << fld.sizes() << ' '
        << flatOutput(fld);

    return Info;
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Basic addressing checks
    #if 0
    {
        ijkAddressing addr1(3, 4, 5);

        Info<< "addressing: " << addr1.sizes() << nl;
        Info<< "index of (2,2,2) " << addr1.index(2,2,2) << nl;

        for (const label idx : labelRange(addr1.size()))
        {
            Info<< "index of " << idx << " => " << addr1.index(idx) << nl;
        }

        for (label k=0; k < addr1.sizes().z(); ++k)
        {
            for (label j=0; j < addr1.sizes().y(); ++j)
            {
                for (label i=0; i < addr1.sizes().x(); ++i)
                {
                    labelVector ijk(i,j,k);

                    Info<< "index of " << addr1.index(ijk)
                        << " <= " << ijk << nl;
                }
            }
        }
    }
    #endif


    // Create with inconsistent sizes
    IjkField<label> field0({1, 2, 3}, identity(10));

    IjkField<scalar> field1({3, 4, 5});
    IjkField<scalar> field2({2, 3, 3});

    forAll(field1, i)
    {
        field1[i] = -i;
    }
    forAll(field2, i)
    {
        field2[i] = i;
    }

    Info<< "ijk field "; print(field1) << nl;
    Info<< "ijk field "; print(field2) << nl;

    field1.resize(field2.sizes());

    Info<< "resized "; print(field1) << nl;

    field1 *= 2;

    Info<< "Multiply: "; print(field1) << nl;

    field1.resize({1, 2, 3});

    Info<< "Resize - shrink: "; print(field1) << nl;

    field1.resize({2, 3, 2});

    Info<< "Resize - grow: "; print(field1) << nl;

    field1.resize({3, 2, 2});

    Info<< "Resize - repartition: "; print(field1) << nl;

    field1 = field2;

    Info<< "Copied: "; print(field1) << nl;

    field1 = 3.14159;

    Info<< "Assigned: "; print(field1) << nl;

    field1 += 3.14159;

    Info<< "+= operator: "; print(field1) << nl;

    field1 /= 1.2;

    Info<< "/= operator: "; print(field1) << nl;

    IjkField<scalar> field3(std::move(field2));

    Info<< "Move construct: "; print(field2) << nl;
    print(field3) << nl;

    // Field operations are still limited, but we can bypass things too

    {
        Field<scalar>& tmpField = field1;
        tmpField = sqr(tmpField);

        Info<< "squared (workaround): "; print(field1) << nl;
    }


    Info<< nl
        << "Before transfer: addr:" << name(field1.data())
        << " size:" << field1.size() << nl;

    Field<scalar> sfield1(std::move(field1));
    field1.clear();

    Info<< "After transfer to regular field" << nl
        << "    source:" << name(field1.data()) << nl
        << "    target:" << name(sfield1.data()) << nl
        << "Values"
        << "    source:";
    print(field1) << nl;

    Info<< "    target:" << flatOutput(sfield1) << nl;


    Info << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //

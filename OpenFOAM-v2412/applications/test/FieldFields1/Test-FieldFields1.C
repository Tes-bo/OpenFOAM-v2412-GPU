/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Application
    Test-FieldFields1

\*---------------------------------------------------------------------------*/

#include "symmTensorField.H"
#include "tensorField.H"
#include "FieldFields.H"
#include "Random.H"

using namespace Foam;


template<class Cmpt>
void printFieldField(const FieldField<Field, Cmpt>& ff)
{
    forAll(ff, i)
    {
        Info<< i << ": " << flatOutput(ff[i]) << nl;
    }
    Info<< nl;
}


template<class Type>
tmp<Field<Type>> randomField(Random& rnd, label dim)
{
    auto tfld = tmp<Field<Type>>::New(dim);
    auto& fld = tfld.ref();

    for (Type& val : fld)
    {
        for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; ++cmpt)
        {
            setComponent(val, cmpt) = rnd.position<label>(0, 100);
        }
    }

    return tfld;
}


template<class Type>
tmp<Field<Type>> randomField(label dim)
{
    Random rnd;
    return randomField<Type>(rnd, dim);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // scalarField
    {
        Info<< nl << "scalarFieldField" << nl;

        Random rnd;

        FieldField<Field, scalar> sff1(6);
        forAll(sff1, i)
        {
            sff1.set(i, randomField<scalar>(rnd, 8));
        }

        printFieldField(sff1);

        Info<< nl << "indexing:" << nl;

        {
            labelPair index;
            const label range1 = sff1.size()-1;
            const label range2 = sff1[0].size()-1;

            for (label iter = 0; iter < 10; ++iter)
            {
                index.first() = rnd.position<label>(0, range1);
                index.second() = rnd.position<label>(0, range2);

                Info<< index << " => " << sff1[index] << nl;
            }
        }
    }

    Info<< nl << "End\n" << nl;

    return 0;
}


// ************************************************************************* //

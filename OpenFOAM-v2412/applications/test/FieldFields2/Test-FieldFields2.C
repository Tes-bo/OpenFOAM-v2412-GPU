/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Test-FieldFields2

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
        Info<< ' ' << flatOutput(ff[i]);
    }
    Info<< nl;
}


template<class Cmpt, unsigned N>
void printComponents(const FixedList<FieldField<Field, Cmpt>, N>& cmpts)
{
    forAll(cmpts, cmpti)
    {
        Info<< cmpti << ':';
        printFieldField(cmpts[cmpti]);
    }
}


template<class Cmpt, unsigned N>
void allocComponents
(
    FixedList<FieldField<Field, Cmpt>, N>& cmpts,
    label dim = 1
)
{
    forAll(cmpts, i)
    {
        cmpts[i] = FieldField<Field, Cmpt>(1);
        cmpts[i].set(0, new Field<Cmpt>(dim, Zero));
    }
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
    const tensor t19(1, 2, 3, 4, 5, 6, 7, 8, 9);

    // tensorField
    {
        Info<< nl << "tensorFieldField" << nl;

        FieldField<Field, tensor> tff1(1);
        tff1.set(0, new tensorField(one{}, t19));

        FixedList<FieldField<Field, scalar>, 9> cmpts;
        allocComponents(cmpts, 1);

        FixedList<FieldField<Field, vector>, 3> slice;
        allocComponents(slice, 1);

        Info<< "=>";
        printFieldField(tff1);

        unzip
        (
            tff1,
            cmpts[0], cmpts[1], cmpts[2],
            cmpts[3], cmpts[4], cmpts[5],
            cmpts[6], cmpts[7], cmpts[8]
        );

        Info<< "components" << nl;
        printComponents(cmpts);

        // Transposed
        zip
        (
            tff1,
            cmpts[0], cmpts[3], cmpts[6],
            cmpts[1], cmpts[4], cmpts[7],
            cmpts[2], cmpts[5], cmpts[8]
        );

        Info<< "rezip (transposed): "
            << " => " << tff1 << nl;


        // Col / rows

        Info<< nl;

        unzipRows(tff1, slice[0], slice[1], slice[2]);
        Info<< "rows" << nl;
        printComponents(slice);

        unzipCols(tff1, slice[0], slice[1], slice[2]);
        Info<< "cols" << nl;
        printComponents(slice);

        // forAll(slice, i)
        // {
        //     unzipRow(tff1, vector::components(i), slice[i]);
        //     Info<< "row " << i << ": " << slice[i] << nl;
        //
        //     unzipCol(tff1, vector::components(i), slice[i]);
        //     Info<< "col " << i << ": " << slice[i] << nl;
        // }

        // Treat columns like rows
        zipCols(tff1, slice[0], slice[1], slice[2]);

        Info<<"rezip (re-transposed): "
            <<" => " << tff1 << nl;

        zipRows(tff1, slice[0], slice[1], slice[2]);
        Info<<"rezip (regular-transposed): "
            <<" => " << tff1 << nl;

        unzipDiag(tff1, slice[0]);
        Info<< "diag" << nl;
        printFieldField(slice[0]);
    }


    // symmTensorField
    {
        Info<< nl << "symmTensorField" << nl;

        FieldField<Field, symmTensor> sf1(1);
        sf1.set(0, new symmTensorField(one{}, symmTensor(1, 2, 3, 4, 5, 6)));

        FixedList<FieldField<Field, scalar>, 6> cmpts;
        allocComponents(cmpts, 1);

        FixedList<FieldField<Field, vector>, 1> slice;
        allocComponents(slice, 1);

        Info<<" =>";
        printFieldField(sf1);

        Info<< nl;
        for (direction cmpt = 0; cmpt < pTraits<vector>::nComponents; ++cmpt)
        {
            unzipRow(sf1, cmpt, slice[0]);

            Info<< "row " << label(cmpt) << ": ";
            printFieldField(slice[0]);
        }

        unzipDiag(sf1, slice[0]);

        Info<< "diag: ";
        printFieldField(slice[0]);

        Info<< nl;

        unzip
        (
            sf1,
            cmpts[0], cmpts[1], cmpts[2],
            cmpts[3], cmpts[4],
            cmpts[5]
        );

        Info<< "components" << nl;
        printComponents(cmpts);

        // Transposed
        zip
        (
            sf1,
            cmpts[5], cmpts[1], cmpts[2],
            cmpts[3], cmpts[4],
            cmpts[1]
        );

        Info<< "rezip (swapped diag): "
            << " => " << sf1 << nl;
    }

    // sphericalTensorField
    {
        Info<< nl << "sphericalTensorField" << nl;

        FieldField<Field, sphericalTensor> sf1(1);
        sf1.set(0, new sphericalTensorField(one{}, sphericalTensor(4)));

        FixedList<FieldField<Field, scalar>, 1> cmpts;
        allocComponents(cmpts, 1);

        Info<<" =>";
        printFieldField(sf1);

        unzip(sf1, cmpts[0]);

        Info<< "components" << nl;
        printComponents(cmpts);
    }

    // vectorField
    {
        Info<< nl << "vectorField" << nl;

        Random rnd;
        FieldField<Field, vector> vf1(1);
        vf1.set(0, randomField<vector>(rnd, 4));

        FixedList<FieldField<Field, scalar>, 3> cmpts;
        allocComponents(cmpts, 4);

        Info<< nl
            << " => " << vf1 << nl;

        unzip(vf1, cmpts[0], cmpts[1], cmpts[2]);

        Info<<"unzip:" << nl
            << "x : " << cmpts[0] << nl
            << "y : " << cmpts[1] << nl
            << "z : " << cmpts[2] << nl
            << nl;

        // Transposed
        zip(vf1, cmpts[2], cmpts[0], cmpts[1]);

        Info<<"rezip (rotated): "
            <<" => " << vf1 << nl;
    }

    Info<< nl << "End\n" << nl;

    return 0;
}


// ************************************************************************* //

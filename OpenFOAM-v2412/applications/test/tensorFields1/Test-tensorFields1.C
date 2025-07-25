/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Test-tensorFields1

\*---------------------------------------------------------------------------*/

#include "tensorField.H"
#include "Random.H"

using namespace Foam;

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
    // tensorField
    {
        Info<< nl << "tensorField" << nl;
        tensorField tf1(2, tensor(1, 2, 3, 4, 5, 6, 7, 8, 9));
        tf1.last() = tf1.last().T();

        FixedList<scalarField, 9> cmpts(scalarField(tf1.size()));

        Info<<" => " << tf1 << nl;

        Info<< nl
            << "row 0: " << unzipRow(tf1, vector::X) << nl
            << "row 1: " << unzipRow(tf1, vector::Y) << nl
            << "row 2: " << unzipRow(tf1, vector::Z) << nl;

        Info<< nl
            << "col 0: " << unzipCol(tf1, vector::X) << nl
            << "col 1: " << unzipCol(tf1, vector::Y) << nl
            << "col 2: " << unzipCol(tf1, vector::Z) << nl;

        Info<< nl
            << "diag: " << unzipDiag(tf1) << nl;

        unzip
        (
            tf1,
            cmpts[0], cmpts[1], cmpts[2],
            cmpts[3], cmpts[4], cmpts[5],
            cmpts[6], cmpts[7], cmpts[8]
        );

        Info<<"unzip:" << nl
            << "xx : " << cmpts[0] << nl
            << "xy : " << cmpts[1] << nl
            << "xz : " << cmpts[2] << nl
            << "yx : " << cmpts[3] << nl
            << "yy : " << cmpts[4] << nl
            << "yz : " << cmpts[5] << nl
            << "zx : " << cmpts[6] << nl
            << "zy : " << cmpts[7] << nl
            << "zz : " << cmpts[8] << nl
            << nl;

        // Transposed
        zip
        (
            tf1,
            cmpts[0], cmpts[3], cmpts[6],
            cmpts[1], cmpts[4], cmpts[7],
            cmpts[2], cmpts[5], cmpts[8]
        );

        Info<<"rezip (transposed): "
            <<" => " << tf1 << nl;

        // Col / rows
        FixedList<vectorField, 3> slice(vectorField(tf1.size()));

        unzipRows(tf1, slice[0], slice[1], slice[2]);

        Info<< nl
            << "rows" << nl
            << " 0: " << slice[0] << nl
            << " 1: " << slice[1] << nl
            << " 2: " << slice[2] << nl;

        unzipCols(tf1, slice[0], slice[1], slice[2]);
        Info<< nl
            << "cols" << nl
            << " 0: " << slice[0] << nl
            << " 1: " << slice[1] << nl
            << " 2: " << slice[2] << nl;


        // Treat columns like rows
        zipCols(tf1, slice[0], slice[1], slice[2]);

        Info<<"rezip (re-transposed): "
            <<" => " << tf1 << nl;

        zipRows(tf1, slice[0], slice[1], slice[2]);
        Info<<"rezip (regular-transposed): "
            <<" => " << tf1 << nl;
    }


    // symmTensorField
    {
        Info<< nl << "symmTensorField" << nl;

        tensorField tf1(2, tensor(1, 2, 3, 4, 5, 6, 7, 8, 9));
        tf1.last() = tf1.last().T();
        scalarField f1(2, one{});

        symmTensorField sf1(2, symmTensor::one);
        symmTensorField sf2(2, symmTensor::one);

        Info<< (tf1 & sf2) << nl;

        Info<< f1*sf1 << " " << sf1*3 << nl;

        Info<< ((sf1 + sf2) & (sf1 + sf2)) << nl;

        vectorField vf1(1, vector::one);
        Info<< sqr(vf1) << nl;
        Info<< pow<vector, 2>(vf1) << nl;

        sf1 = symm(tf1);

        Info<< sf1 << nl;
        Info<<" => " << tf1 << nl;

        sf1 = symmTensor(1, 2, 3, 4, 5, 6);
        Info<< nl << "reset: " << sf1 << nl;

        FixedList<scalarField, 6> cmpts(scalarField(sf1.size()));

        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            Info<< "col/row " << label(cmpt) << ": "
                << sf1[0].col(cmpt) << ' '
                << sf1[0].row(cmpt) << nl;
        }

        Info<< nl
            << "diag: " << unzipDiag(sf1) << nl
            << "row 0: " << unzipRow(sf1, vector::X) << nl
            << "row 1: " << unzipCol(sf1, vector::Y) << nl  // same as row
            << "row 2: " << unzipRow(sf1, vector::Z) << nl;


        unzip
        (
            sf1,
            cmpts[0], cmpts[1], cmpts[2],
            cmpts[3], cmpts[4],
            cmpts[5]
        );

        Info<< "unzip:" << nl
            << "xx : " << cmpts[0] << nl
            << "xy : " << cmpts[1] << nl
            << "xz : " << cmpts[2] << nl
            << "yy : " << cmpts[3] << nl
            << "yz : " << cmpts[4] << nl
            << "zz : " << cmpts[5] << nl
            << nl;

        // Transposed diagonal
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

        scalarField f1(2, one{});
        sphericalTensorField sf1(2, sphericalTensor(1));
        sphericalTensorField sf2(2, sphericalTensor(2));
        tensorField tf1(2, tensor::one);

        Info<< (tf1 & sf2) << nl;

        Info<< f1*sf1 << " " << sf1*3 << nl;
        Info<< ((sf1 + sf2) & (sf1 + sf2)) << nl;

        sf1[0] = sphericalTensor(2);
        sf1[1] = sphericalTensor(1);
        scalarField cmpts(sf1.size());

        unzip(sf1, cmpts);

        Info<< nl
            <<" => " << sf1 << nl;

        Info<<"unzip:" << nl
            << cmpts << nl
            << nl;
    }


    // vectorField
    {
        Info<< nl << "vectorField" << nl;

        vectorField vf1(randomField<vector>(4));
        FixedList<scalarField, 3> cmpts(scalarField(vf1.size()));

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

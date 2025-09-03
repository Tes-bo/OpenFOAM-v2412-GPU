/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "MatrixTools.H"

// * * * * * * * * * * * * * * * Global Functions * * * * * * * * * * * * * * //

template<class Form1, class Form2, class Type>
bool Foam::MatrixTools::equal
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B,
    const bool verbose,
    const label maxDiffs,
    const scalar relTol,
    const scalar absTol
)
{
    const label len = A.size();

    if (len != B.size())
    {
        if (verbose)
        {
            Info<< "Matrices have different sizes: "
                << len << " vs " << B.size() << nl;
        }
        return false;
    }

    auto iter1 = A.cbegin();
    auto iter2 = B.cbegin();

    label nDiffs = 0;

    for (label i = 0; i < len; ++i)
    {
        if ((absTol + relTol*mag(*iter2)) < Foam::mag(*iter1 - *iter2))
        {
            ++nDiffs;

            if (verbose)
            {
                Info<< "Matrix element " << i
                    << " differs beyond tolerance: "
                    << *iter1 << " vs " << *iter2 << nl;
            }
            if (maxDiffs && maxDiffs < nDiffs)
            {
                Info<< "More than " << maxDiffs << " elements differ."
                    << " Skipping further comparisons." << nl;
                return false;
            }
        }

        ++iter1;
        ++iter2;
    }

    if (verbose)
    {
        if (nDiffs)
        {
            Info<< "Some elements differed" << nl;
        }
        else
        {
            Info<< "All elements equal within the tolerances" << nl;
        }
    }

    return !nDiffs;
}


template<class Container>
Foam::Ostream& Foam::MatrixTools::printMatrix
(
    Ostream& os,
    const Container& mat
)
{
    os  << mat.m() << ' ' << mat.n();

    if (mat.m() == 1)
    {
        // row-vector
        os  << " (";
        for (label j = 0; j < mat.n(); ++j)
        {
            if (j) os  << ' ';
            os  << mat(0,j);
        }
        os  << ')' << nl;
    }
    else if (mat.n() == 1)
    {
        // col-vector

        os  << " (";
        for (label i = 0; i < mat.m(); ++i)
        {
            if (i) os  << ' ';
            os  << mat(i,0);
        }
        os  << ')' << nl;
    }
    else
    {
        // Regular

        os  << nl << '(' << nl;

        for (label i = 0; i < mat.m(); ++i)
        {
            os  << '(';
            for (label j = 0; j < mat.n(); ++j)
            {
                if (j) os  << ' ';
                os  << mat(i,j);
            }
            os  << ')' << nl;
        }
        os  << ')' << nl;
    }

    return os;
}


// ************************************************************************* //

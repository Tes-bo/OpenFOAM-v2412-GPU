/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-invTensor

Description
    Tests for regular and corner cases of tensor inversion.

\*---------------------------------------------------------------------------*/

#include "tensor.H"
#include "symmTensor.H"
#include "transform.H"
#include "unitConversion.H"
#include "Random.H"
#include "scalar.H"
#include "complex.H"
#include "sigFpe.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TensorType>
void test_invert(const TensorType& tt)
{
    Info<< pTraits<TensorType>::typeName << nl
        << "ten : " << tt << nl;

    // Non-failsafe. try/catch does not work here
    if (tt.det() < ROOTVSMALL)
    {
        Info<< "inv : " << "nan/inf" << nl;
    }
    else
    {
        Info<< "inv : " << tt.inv() << nl;
    }

    // Failsafe
    Info<< "inv : " << tt.safeInv() << nl;

    Info<< nl;
}


// * * * * * * * * * * * * * * * Main Program  * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Run with FOAM_SIGFPE=true to ensure we see divide by zero
    Foam::sigFpe::set(true);

    {
        symmTensor st(1e-3, 0, 0, 1e-3, 0, 1e-6);
        test_invert(st);
    }

    {
        symmTensor st(1e-3, 0, 0, 1e-3, 0, 1e-30);
        test_invert(st);
    }

    {
        symmTensor st(1e-3, 0, 0, 1e-3, 0, 0);
        test_invert(st);
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "fft.H"
#include <fftw3.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::fft::fftRenumberRecurse
(
    List<complex>& data,
    List<complex>& renumData,
    const UList<int>& nn,
    label nnprod,
    label ii,
    label l1,
    label l2
)
{
    if (ii == nn.size())
    {
        // We've worked out the renumbering scheme. Now copy
        // the components across

        data[l1] = renumData[l2];
    }
    else
    {
        // Do another level of folding. First work out the
        // multiplicative value of the index

        nnprod /= nn[ii];
        label i_1(0);

        for (label i=0; i<nn[ii]; i++)
        {
            // Now evaluate the indices (both from array 1 and to
            // array 2). These get multiplied by nnprod to (cumulatively)
            // find the real position in the list corresponding to
            // this set of indices

            if (i < nn[ii]/2)
            {
                i_1 = i + nn[ii]/2;
            }
            else
            {
                i_1 = i - nn[ii]/2;
            }


            // Go to the next level of recursion

            fftRenumberRecurse
            (
                data,
                renumData,
                nn,
                nnprod,
                ii+1,
                l1+i*nnprod,
                l2+i_1*nnprod
            );
        }
    }
}


void Foam::fft::fftRenumber(List<complex>& data, const UList<int>& nn)
{
    List<complex> renumData(data);

    label nnprod(1);
    forAll(nn, i)
    {
        nnprod *= nn[i];
    }

    label ii(0), l1(0), l2(0);

    fftRenumberRecurse
    (
        data,
        renumData,
        nn,
        nnprod,
        ii,
        l1,
        l2
    );
}


Foam::tmp<Foam::complexField>
Foam::fft::realTransform1D(const scalarField& field)
{
    const label n = field.size();
    const label nBy2 = n/2;

    // Copy of input field for use by fftw
    // - require non-const access to input and output
    // - use double to avoid additional libfftwf for single-precision

    List<double> in(n);
    List<double> out(n);

    for (label i=0; i < n; ++i)
    {
        in[i] = field[i];
    }

    // Using real to half-complex fftw 'kind'
    fftw_plan plan = fftw_plan_r2r_1d
    (
        n,
        in.data(),
        out.data(),
        FFTW_R2HC,
        FFTW_ESTIMATE
    );

    fftw_execute(plan);

    // field[0] = DC component
    auto tresult = tmp<complexField>::New(nBy2 + 1);
    auto& result = tresult.ref();

    result[0].real(out[0]);
    result[nBy2].real(out[nBy2]);
    for (label i = 1; i < nBy2; ++i)
    {
        result[i].real(out[i]);
        result[i].imag(out[n - i]);
    }

    fftw_destroy_plan(plan);

    return tresult;
}


Foam::tmp<Foam::complexField> Foam::fft::realTransform1D
(
    const tmp<scalarField>& tfield
)
{
    tmp<complexField> tresult = realTransform1D(tfield());
    tfield.clear();
    return tresult;
}


void Foam::fft::transform
(
    complexField& field,
    const UList<int>& nn,
    transformDirection dir
)
{
    // Copy field into fftw containers
    const label N = field.size();

    fftw_complex* inPtr =
        static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N));
    fftw_complex* outPtr =
        static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*N));

    // If reverse transform : renumber before transform
    if (dir == REVERSE_TRANSFORM)
    {
        fftRenumber(field, nn);
    }

    forAll(field, i)
    {
        inPtr[i][0] = field[i].real();
        inPtr[i][1] = field[i].imag();
    }

    // Create the plan
    // FFTW_FORWARD = -1
    // FFTW_BACKWARD = 1

    // 1-D plan
    // fftw_plan plan =
    //    fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Generic 1..3-D plan
    const label rank = nn.size();
    fftw_plan plan =
        fftw_plan_dft(rank, nn.begin(), inPtr, outPtr, dir, FFTW_ESTIMATE);

    // Compute the FFT
    fftw_execute(plan);

    forAll(field, i)
    {
        field[i].real(outPtr[i][0]);
        field[i].imag(outPtr[i][1]);
    }

    fftw_destroy_plan(plan);

    fftw_free(inPtr);
    fftw_free(outPtr);

    // If forward transform : renumber after transform
    if (dir == FORWARD_TRANSFORM)
    {
        fftRenumber(field, nn);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::complexField> Foam::fft::forwardTransform
(
    const tmp<complexField>& tfield,
    const UList<int>& nn
)
{
    auto tresult = tmp<complexField>::New(tfield);

    transform(tresult.ref(), nn, FORWARD_TRANSFORM);

    tfield.clear();

    return tresult;
}


Foam::tmp<Foam::complexField> Foam::fft::reverseTransform
(
    const tmp<complexField>& tfield,
    const UList<int>& nn
)
{
    auto tresult = tmp<complexField>::New(tfield);

    transform(tresult.ref(), nn, REVERSE_TRANSFORM);

    tfield.clear();

    return tresult;
}


Foam::tmp<Foam::complexVectorField> Foam::fft::forwardTransform
(
    const tmp<complexVectorField>& tfield,
    const UList<int>& nn
)
{
    auto tresult = tmp<complexVectorField>::New(tfield().size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        tresult.ref().replace
        (
            cmpt,
            forwardTransform(tfield().component(cmpt), nn)
        );
    }

    tfield.clear();

    return tresult;
}


Foam::tmp<Foam::complexVectorField> Foam::fft::reverseTransform
(
    const tmp<complexVectorField>& tfield,
    const UList<int>& nn
)
{
    auto tresult = tmp<complexVectorField>::New(tfield().size());

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        tresult.ref().replace
        (
            cmpt,
            reverseTransform(tfield().component(cmpt), nn)
        );
    }

    tfield.clear();

    return tresult;
}


// ************************************************************************* //

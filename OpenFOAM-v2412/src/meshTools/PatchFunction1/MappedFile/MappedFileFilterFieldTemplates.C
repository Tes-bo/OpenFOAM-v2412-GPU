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

\*---------------------------------------------------------------------------*/

#include "indexedOctree.H"
#include "treeDataPoint.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::PatchFunction1Types::FilterField::evaluate
(
    const tmp<Field<Type>>& tinput,
    const label nSweeps
) const
{
    if (nSweeps < 1 || !tinput.good())
    {
        // Nothing to do
        return tinput;
    }

    label nPoints = tinput().size();

    const label nAddr = addressing_.size();

    if (!nPoints || !nAddr)
    {
        // No input or using an identity mapping
        return tinput;
    }

    // Output
    auto toutput = tmp<Field<Type>>::New(nPoints);

    if (nAddr < nPoints)
    {
        WarningInFunction
            << "Addressing/weights shorter than input field"
            << endl;

        // Restrict addressing space within loop
        nPoints = nAddr;

        // Fill trailing portion with the input values
        toutput.ref().slice(nAddr) = tinput().slice(nAddr);
    }

    // Intermediate buffer
    tmp<Field<Type>> tbuffer;

    if (nSweeps == 1)
    {
        // Simple reference ie enough
        tbuffer.cref(tinput);
    }
    else
    {
        // Need buffer for swap - copy or steal contents
        tbuffer.reset(tinput.ptr());
    }
    tinput.clear();

    // If there are any 'senseless' sweeps
    // (ie, with identity mapping of values), they will have
    // been detected during addressing construction, can ignore here

    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        if (sweep > 0)
        {
            // Reuse previous output for subsequent sweeps
            tbuffer.swap(toutput);
        }

        const auto& input = tbuffer();
        auto& output = toutput.ref();

        #pragma omp parallel for if (nPoints > 1000)
        for (label pointi = 0; pointi < nPoints; ++pointi)
        {
            const auto& addr = addressing_[pointi];
            const auto& weight = weights_[pointi];

            auto& interp = output[pointi];

            if (addr.empty())
            {
                // Could happen if the search radius was really small
                interp = input[pointi];
            }
            else
            {
                interp = Zero;

                forAll(addr, i)
                {
                    interp += (weight[i] * input[addr[i]]);
                }
            }
        }
    }

    return toutput;
}


// ************************************************************************* //

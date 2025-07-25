/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "deltaBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class pT>
pT deltaBoundary::makeCellCentres_d
(
    const vectorField& fAreas,
    const vectorField& fCtrs,
    const Field<pT>& fAreas_d,
    const Field<pT>& fCtrs_d
)
{
    // Define type that in an order smaller than pT. Used for volume-related
    // variations
    typedef typename innerProduct<vector, pT>::type vT;

    // First estimate the approximate cell centre as the average of
    // face centres
    vector cEst(Zero);
    vector cellCtrs(Zero);
    scalar cellVols(Zero);
    pT cEst_d(pTraits<pT>::zero);
    pT cellCtrs_d(pTraits<pT>::zero);
    vT cellVols_d(pTraits<vT>::zero);

    forAll(fAreas, facei)
    {
        cEst += fCtrs[facei];
        cEst_d += fCtrs_d[facei];
    }

    cEst /= fAreas.size();
    cEst_d /= fAreas.size();

    forAll(fAreas, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol =
            mag(fAreas[facei] & (fCtrs[facei] - cEst));

        vT pyr3Vol_d =
            (fAreas[facei] & (fCtrs[facei] - cEst))
           *(
               ((fCtrs[facei] - cEst) & fAreas_d[facei])
               // Reverse order to get the correct inner product
             + (fAreas[facei] & (fCtrs_d[facei] - cEst_d))
            )/pyr3Vol;

        // Calculate face-pyramid centre
        vector pc = (3.0/4.0)*fCtrs[facei] + (1.0/4.0)*cEst;
        pT pc_d = (3.0/4.0)*fCtrs_d[facei] + (1.0/4.0)*cEst_d;

        // Accumulate volume-weighted face-pyramid centre
        cellCtrs += pyr3Vol*pc;

        // Reverse order to get the correct outer product
        cellCtrs_d += (pc*pyr3Vol_d + pyr3Vol*pc_d);

        // Accumulate face-pyramid volume
        cellVols += pyr3Vol;
        cellVols_d += pyr3Vol_d;
    }

    cellCtrs /= cellVols;
    cellCtrs_d = cellCtrs_d/cellVols - cellCtrs*cellVols_d/cellVols;
    cellVols *= (1.0/3.0);
    cellVols_d *= (1.0/3.0);

    return cellCtrs_d;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "radialActuationDiskSource.H"
#include "volFields.H"
#include "fvMatrix.H"
#include "fvm.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::radialActuationDiskSource::
addRadialActuationDiskAxialInertialResistance
(
    vectorField& Usource,
    const labelList& cells,
    const scalarField& Vcells,
    const RhoFieldType& rho,
    const vectorField& U
)
{
    scalarField Tr(cells.size());

    tensor E(Zero);
    const vector diskDir = this->diskDir();
    E.diag(diskDir);

    const Field<vector> zoneCellCentres(mesh().cellCentres(), cells);
    const Field<scalar> zoneCellVolumes(mesh().cellVolumes(), cells);

    const vector avgCentre = gSum(zoneCellVolumes*zoneCellCentres)/V();
    const scalar maxR = gMax(mag(zoneCellCentres - avgCentre));

    const scalar intCoeffs =
        radialCoeffs_[0]
      + radialCoeffs_[1]*sqr(maxR)/2.0
      + radialCoeffs_[2]*pow4(maxR)/3.0;

    if (mag(intCoeffs) < VSMALL)
    {
        FatalErrorInFunction
            << "Radial distribution coefficients lead to zero polynomial." << nl
            << "radialCoeffs = " << radialCoeffs_
            << exit(FatalError);
    }

    // Compute upstream U and rho, spatial-averaged over monitor-region
    vector Uref(Zero);
    scalar rhoRef = 0.0;
    label szMonitorCells = monitorCells_.size();

    for (const label celli : monitorCells_)
    {
        Uref += U[celli];
        rhoRef = rhoRef + rho[celli];
    }
    reduce(Uref, sumOp<vector>());
    reduce(rhoRef, sumOp<scalar>());
    reduce(szMonitorCells, sumOp<label>());

    if (szMonitorCells == 0)
    {
        FatalErrorInFunction
            << "No cell is available for incoming velocity monitoring."
            << exit(FatalError);
    }

    Uref /= szMonitorCells;
    rhoRef /= szMonitorCells;

    const scalar Ct = sink_*UvsCtPtr_->value(mag(Uref));
    const scalar Cp = sink_*UvsCpPtr_->value(mag(Uref));

    if (Cp <= VSMALL || Ct <= VSMALL)
    {
        FatalErrorInFunction
           << "Cp and Ct must be greater than zero." << nl
           << "Cp = " << Cp << ", Ct = " << Ct
           << exit(FatalError);
    }

    const scalar a = 1.0 - Cp/Ct;
    const scalar T = 2.0*rhoRef*diskArea_*mag(Uref)*a*(1.0 - a);

    forAll(cells, i)
    {
        const scalar r2 = magSqr(mesh().cellCentres()[cells[i]] - avgCentre);

        Tr[i] =
            T
           *(radialCoeffs_[0] + radialCoeffs_[1]*r2 + radialCoeffs_[2]*sqr(r2))
           /intCoeffs;

        Usource[cells[i]] += ((Vcells[cells[i]]/V_)*Tr[i]*E) & Uref;
    }

    if
    (
        mesh_.time().timeOutputValue() >= writeFileStart_
     && mesh_.time().timeOutputValue() <= writeFileEnd_
    )
    {
        Ostream& os = file();
        writeCurrentTime(os);

        os  << Uref << tab << Cp << tab << Ct << tab << a << tab << T << tab
            << endl;
    }

    if (debug)
    {
        Info<< "Source name: " << name() << nl
            << "Average centre: " << avgCentre << nl
            << "Maximum radius: " << maxR << endl;
    }
}


// ************************************************************************* //

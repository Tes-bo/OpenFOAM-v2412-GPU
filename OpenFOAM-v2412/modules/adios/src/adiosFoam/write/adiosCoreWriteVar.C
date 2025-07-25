/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "adiosCoreWrite.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adiosFoam::adiosCoreWrite::putLabelVariable
(
    const fileName& name,
    const label value
)
{
    using cmptType = label;

    const auto global =
        static_cast<adios2::Dims::value_type>(UPstream::nProcs());

    const auto offset =
        static_cast<adios2::Dims::value_type>(UPstream::myProcNo());


    // Define and write
    auto var =
        writeIOPtr_->DefineVariable<cmptType>
        (
            name,
            { global }, // shape = global 1D array of this info
            { offset }, // start (offset) = for this rank
            { 1 },      // single value per rank
            adios2::ConstantDims // constant shape over time
        );

    writeFilePtr_->Put<cmptType>(var, value, adios2::Mode::Sync);
}


// ************************************************************************* //

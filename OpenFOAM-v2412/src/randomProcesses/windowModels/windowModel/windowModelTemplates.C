/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::windowModel::apply
(
    const Field<Type>& fld,
    const label windowI
) const
{
    label nSamples = this->nSamples();

    if (nSamples > fld.size())
    {
        FatalErrorInFunction
            << "Number of samples in sampling window is greater than the "
            << "size of the input field" << nl
            << "    input field size       = " << fld.size() << nl
            << "    window size            = " << nSamples << nl
            << "    requested window index = " << windowI
            << exit(FatalError);
    }


    auto tresult = tmp<Field<Type>>::New(nSamples, Zero);
    auto& result = tresult.ref();

    label nWindow = nWindowsTotal(fld.size());
    if (windowI >= nWindow)
    {
        FatalErrorInFunction
            << "Requested window " << windowI << " outside of range. "
            << "Number of available windows is " << nWindow
            << abort(FatalError);
    }

    label windowOffset = windowI*(nSamples - nOverlapSamples_);

    const scalarField& wf = *this;
    result = wf*SubField<Type>(fld, nSamples, windowOffset);

    return tresult;
}


// ************************************************************************* //

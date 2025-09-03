/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 DLR
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

#include "leastSquareGrad.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::leastSquareGrad<T>::leastSquareGrad
(
    const word& functionName,
    const labelVector& geomDir
)
:
    polyFitter_(functionName,geomDir),
    geomDir_(geomDir),
    nDims_(0)
{
    // Compute number of dimensions
    for (const label dirn : geomDir_)
    {
        if (dirn == 1)
        {
            ++nDims_;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
typename Foam::outerProduct<Foam::vector, T>::type
Foam::leastSquareGrad<T>::grad
(
    const List<vector>& positions,
    const List<T>& listValue
)
{
    typedef typename outerProduct<vector, T>::type GradType;

    List<T> fitData = polyFitter_.fitData
    (
        positions,
        listValue
    );

    if (nDims_ == 3)
    {
        return GradType(fitData[1],fitData[2],fitData[3]);
    }


    label dimCounter = 0;

    GradType ret(Zero);

    forAll(geomDir_,i)
    {
        if (geomDir_[i] == 1)
        {
            ++dimCounter;
            ret[i] = fitData[dimCounter];
        }
    }

    return ret;
}


namespace Foam // needed g++ bug
{
    template<>
    tensor leastSquareGrad<vector>::grad
    (
        const List<vector>& positions,
        const List<vector>& listValue
    )
    {
        typedef tensor GradType;

        List<vector> fitData = polyFitter_.fitData
        (
            positions,
            listValue
        );

        if (nDims_ == 3)
        {
            return GradType(fitData[1],fitData[2],fitData[3]);
        }

        label dimCounter = 0;

        GradType ret(Zero);

        forAll(geomDir_,i)
        {
            if (geomDir_[i] == 1)
            {
                ++dimCounter;
                ret.row(i, fitData[dimCounter]);
            }
        }

        return ret;
    }
}


template<class T>
Foam::Map<typename Foam::outerProduct<Foam::vector, T>::type>
Foam::leastSquareGrad<T>::grad
(
    const Map<List<vector>>& positions,
    const Map<List<T>>& listValue
)
{
    typedef typename outerProduct<vector, T>::type GradType;

    Map<GradType> gradMap(positions.capacity());

    forAllConstIters(positions, iter)
    {
        const label key = iter.key();
        const List<vector>& positions = iter.val();

        GradType grad(this->grad(positions, listValue[key]));

        gradMap.insert(key, grad);
    }

    return gradMap;
}


template class Foam::leastSquareGrad<Foam::scalar>;
template class Foam::leastSquareGrad<Foam::vector>;


// ************************************************************************* //

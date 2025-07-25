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

Typedef
    Foam::dimensionedScalarMinMax

Description
    A dimensioned scalarMinMix (MinMax for scalar quantities).

Typedef
    Foam::dimensionedMinMax\<T\>

Description
    A templated type alias for dimensioned\<MinMax\<T\>\>

SourceFiles
    dimensionedMinMaxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dimensionedMinMax_H
#define Foam_dimensionedMinMax_H

#include "dimensionedType.H"
#include "MinMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// General alias
template<class T>
using dimensionedMinMax = dimensioned<MinMax<T>>;

// Common typedefs
typedef dimensioned<scalarMinMax> dimensionedScalarMinMax;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Make a dimensionedMinMax from all components
template<class T>
inline dimensioned<MinMax<T>> makeDimensionedMinMax
(
    const word& name,
    const dimensionSet& dims,
    const T& minVal,
    const T& maxVal
)
{
    return dimensioned<MinMax<T>>(name, dims, MinMax<T>(minVal, maxVal));
}


//- Create a dimensionless "minmax"
template<class T>
inline dimensioned<MinMax<T>> makeDimensionedMinMax
(
    const T& minVal,
    const T& maxVal
)
{
    return
        dimensioned<MinMax<T>>
        (
            "minmax", dimensionSet(), MinMax<T>(minVal, maxVal)
        );
}


//- Combine two dimensioned types into a dimensionedMinMax with specified name
template<class T>
inline dimensioned<MinMax<T>> makeDimensionedMinMax
(
    const word& name,
    const dimensioned<T>& minVal,
    const dimensioned<T>& maxVal
)
{
    // Dimension check when dimensionSet::checking()
    return dimensioned<MinMax<T>>
    (
        name,
        (minVal.dimensions() + maxVal.dimensions()),
        MinMax<T>(minVal.value(), maxVal.value())
    );
}


//- Combine two dimensioned types into a dimensionedMinMax "minmax"
template<class T>
inline dimensioned<MinMax<T>> makeDimensionedMinMax
(
    const dimensioned<T>& minVal,
    const dimensioned<T>& maxVal
)
{
    return makeDimensionedMinMax("minmax", minVal, maxVal);
}


//- Construct from components (name, dimensions, value) with an optional
//- dictionary override that can also \e zip together different sub-entries.
//
//  The dictionary override can specify separate min/max dictionary entries.
//  For example,
//  \verbatim
//      rhoMin          0.3;
//      rhoMax          2.0;
//  \endverbatim
//
//  Construct as following:
//
//  \verbatim
//      makeDimensionedMinMax<scalar>
//      (
//          "rhoMinMax", dimDensity, {Zero, GREAT}, dict, "rhoMin", "rhoMax"
//      );
//  \endverbatim
template<class T>
dimensioned<MinMax<T>> makeDimensionedMinMax
(
    const word& name,
    const dimensionSet& dims,
    const MinMax<T>& values,
    const dictionary& dict,
    const word& minName,
    const word& maxName
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dimensionedMinMaxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

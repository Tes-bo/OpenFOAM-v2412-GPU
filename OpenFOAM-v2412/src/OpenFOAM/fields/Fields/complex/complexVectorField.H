/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
    Foam::complexVectorField

Description
    Specialisation of Field\<T\> for complexVector.

SourceFiles
    complexVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_complexVectorField_H
#define Foam_complexVectorField_H

#include "complex.H"
#include "complexVector.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Field<complexVector> complexVectorField;


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Zip together complex vector field from real/imag components
void zip
(
    complexVectorField& result,
    const UList<vector>& realValues,
    const UList<vector>& imagValues
);


//- Zip together complex vector field from real components
//- and constant imag component
void zip
(
    complexVectorField& result,
    const UList<vector>& realValues,
    const vector& imagValue
);


//- Zip together complex vector field from constant real component
//- and imag components
void zip
(
    complexVectorField& result,
    const vector& realValue,
    const UList<vector>& imagValues
);


//- Create complexVector field by zipping two lists of real/imag values
complexVectorField ComplexField
(
    const UList<vector>& realValues,
    const UList<vector>& imagValues
);

//- Create complexVector field by zipping a list of real components
//- and a constant imag component
complexVectorField ComplexField
(
    const UList<vector>& realValues,
    const vector& imagValue
);

//- Create complexVector field by zipping a constant real component
//- and a list of imag components
complexVectorField ComplexField
(
    const vector& realValue,
    const UList<vector>& imagValues
);


//- Extract real component
vectorField Re(const UList<complexVector>& cmplx);

//- Extract imag component
vectorField Im(const UList<complexVector>& cmplx);

//- Sum real and imag components
vectorField ReImSum(const UList<complexVector>& cmplx);


//- Create complexVector field from a list of real (using imag == 0)
inline complexVectorField ReComplexField(const UList<vector>& realValues)
{
    return ComplexField(realValues, vector::zero);
}

//- Create complexVector field from a list of imag (using real == 0)
inline complexVectorField ImComplexField(const UList<vector>& imagValues)
{
    return ComplexField(vector::zero, imagValues);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Cross-product
complexVectorField operator^
(
    const UList<vector>& vec,
    const UList<complexVector>& cmplx
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

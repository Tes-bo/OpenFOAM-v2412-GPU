/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

InNamespace
    Foam::ensightOutput

Description
    A collection of functions for writing volField content in ensight format.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightOutput_volField_H
#define Foam_ensightOutput_volField_H

#include "ensightOutput.H"
#include "ensightFaces.H"
#include "ensightCells.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ensightMesh;

namespace ensightOutput
{

/*---------------------------------------------------------------------------*\
                            Namespace ensightOutput
\*---------------------------------------------------------------------------*/

//- Write volume field component-wise
template<class Type>
bool writeVolField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, fvPatchField, volMesh>& vf,

    //! The addressing (element-wise) into the field
    const ensightMesh& ensMesh
);

//- Write volume field component-wise
template<class Type>
bool writeVolField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, fvPatchField, volMesh>& vf,

    //! The addressing (element-wise) into the field
    const ensightMesh& ensMesh
)
{
    ensightOutput::floatBufferType scratch;
    return ensightOutput::writeVolField<Type>(scratch, os, vf, ensMesh);
}

//- Write volume field component-wise, optionally forcing interpolation
template<class Type>
bool writeVolField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, fvPatchField, volMesh>& vf,

    //! The addressing (element-wise) into the field
    const ensightMesh& ensMesh,

    //! Interpolate values to points
    const bool nodeValues
);

//- Write volume field component-wise, optionally forcing interpolation
template<class Type>
bool writeVolField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, fvPatchField, volMesh>& vf,

    //! The addressing (element-wise) into the field
    const ensightMesh& ensMesh,

    //! Interpolate values to points
    const bool nodeValues
)
{
    ensightOutput::floatBufferType scratch;
    return ensightOutput::writeVolField<Type>
    (
        scratch, os, vf, ensMesh, nodeValues
    );
}

//- Write point field component-wise
template<class Type>
bool writePointField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, pointPatchField, pointMesh>& pf,

    //! The addressing (element-wise) into the field
    const ensightMesh& ensMesh
);

//- Write point field component-wise
template<class Type>
bool writePointField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, pointPatchField, pointMesh>& pf,

    //! The addressing (element-wise) into the field
    const ensightMesh& ensMesh
)
{
    ensightOutput::floatBufferType scratch;
    return ensightOutput::writePointField<Type>(scratch, os, pf, ensMesh);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputVolFieldTemplates.C"
#endif

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    A collection of functions for writing areaField content in ensight format.

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightOutput_areaField_H
#define Foam_ensightOutput_areaField_H

#include "ensightOutput.H"
#include "ensightFaces.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class ensightFaMesh;

namespace ensightOutput
{

/*---------------------------------------------------------------------------*\
                            Namespace ensightOutput
\*---------------------------------------------------------------------------*/

//- Write finite-area field component-wise
template<class Type>
bool writeAreaField
(
    //! Component scratch buffer
    ensightOutput::floatBufferType& scratch,

    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, faPatchField, areaMesh>& fld,

    //! The addressing into the field
    const ensightFaMesh& ensMesh
);


//- Write finite-area field component-wise
template<class Type>
bool writeAreaField
(
    //! Output file (must be valid on master)
    ensightFile& os,

    //! The field content to be written
    const GeometricField<Type, faPatchField, areaMesh>& fld,

    //! The addressing into the field
    const ensightFaMesh& ensMesh
)
{
    ensightOutput::floatBufferType scratch;
    return ensightOutput::writeAreaField(scratch, os, fld, ensMesh);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensightOutput
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightOutputAreaFieldTemplates.C"
#endif

#endif

// ************************************************************************* //

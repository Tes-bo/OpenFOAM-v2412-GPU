/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

Class
    Foam::surfaceReader

Description
    Abstract base class for surface readers with fields.

    Some readers support different input options, these are typically
    specified as 'readOptions' in the containing dictionary.

    \verbatim
    readOptions
    {
        default
        {
            verbose     false;
        }

        ensight
        {
            masterOnly  false;
        }
    }
    \endverbatim

SourceFiles
    surfaceReader.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceReader_H
#define Foam_surfaceReader_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "MeshedSurfaces.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceReader Declaration
\*---------------------------------------------------------------------------*/

class surfaceReader
{
protected:

    //- File name
    fileName fileName_;


public:

    //- Runtime type information
    TypeName("surfaceReader");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceReader,
            fileName,
            (
                const fileName& fName,
                const dictionary& options
            ),
            (fName, options)
        );


    // Helpers

        //- Same as fileFormats::getFormatOptions
        static dictionary formatOptions
        (
            const dictionary& dict,
            const word& formatName,
            const word& entryName = "formatOptions"
        );


    // Selectors

        //- Return a reference to the selected surfaceReader
        static autoPtr<surfaceReader> New
        (
            const word& readType,
            const fileName& fName,
            const dictionary& options = dictionary()
        );


    // Constructors

        //- Construct from fileName
        explicit surfaceReader(const fileName& fName);

        //- Construct from fileName and specified options
        surfaceReader(const fileName& fName, const dictionary& options);


    //- Destructor
    virtual ~surfaceReader() = default;


    // Member Functions

        //- Return a reference to the surface geometry
        virtual const meshedSurface& geometry(const label timeIndex) = 0;

        //- Return a list of the available times
        virtual instantList times() const = 0;

        //- Return a list of the available fields at a given time
        virtual wordList fieldNames(const label timeIndex) const = 0;

        //- Return a scalar field at a given time
        virtual tmp<Field<scalar>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const scalar& refValue = pTraits<scalar>::zero
        ) const = 0;

        //- Return a vector field at a given time
        virtual tmp<Field<vector>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const vector& refValue = pTraits<vector>::zero
        ) const = 0;

        //- Return a sphericalTensor field at a given time
        virtual tmp<Field<sphericalTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const sphericalTensor& refValue = pTraits<sphericalTensor>::zero
        ) const = 0;

        //- Return a symmTensor field at a given time
        virtual tmp<Field<symmTensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const symmTensor& refValue = pTraits<symmTensor>::zero
        ) const = 0;

        //- Return a tensor field at a given time
        virtual tmp<Field<tensor>> field
        (
            const label timeIndex,
            const label fieldIndex,
            const tensor& refValue = pTraits<tensor>::zero
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::foamToVtkReportFields

Description
    Collection of simple static methods for reporting field names
    by category, which is used by foamToVTK.

\*---------------------------------------------------------------------------*/

#ifndef foamToVtkReportFields_H
#define foamToVtkReportFields_H

#include "UPtrList.H"
#include "Ostream.H"
#include "areaFields.H"
#include "volFields.H"
#include "pointFields.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class foamToVtkReportFields Declaration
\*---------------------------------------------------------------------------*/

struct foamToVtkReportFields
{
    template<class GeoField>
    static void print
    (
        const char* msg,
        Ostream& os,
        const UPtrList<const GeoField>& flds
    )
    {
        if (flds.size())
        {
            os  << msg;
            for (const GeoField& fld : flds)
            {
                os  << ' ' << fld.name();
            }
            os  << endl;
        }
    }


    static void print
    (
        const char* msg,
        Ostream& os,
        const wordList& fieldNames
    )
    {
        if (fieldNames.size())
        {
            os  << msg;
            for (const word& fieldName : fieldNames)
            {
                os  << ' ' << fieldName;
            }
            os  << endl;
        }
    }


    template<class FieldType>
    static void print
    (
        const char* msg,
        Ostream& os,
        const IOobjectList& objects
    )
    {
        print(msg, os, objects.sortedNames<FieldType>());
    }


    //- Supported volume field types
    static void volume(Ostream& os, const IOobjectList& objects)
    {
        print<volScalarField>
        (
            "    volScalar    :", os, objects
        );
        print<volVectorField>
        (
            "    volVector    :", os, objects
        );
        print<volSphericalTensorField>
        (
            "    volSphTensor :", os, objects
        );
        print<volSymmTensorField>
        (
            "    volSymTensor :", os, objects
        );
        print<volTensorField>
        (
            "    volTensor    :", os, objects
        );
    }


    //- Supported dimensioned field types
    static void internal(Ostream& os, const IOobjectList& objects)
    {
        print<volScalarField::Internal>
        (
            "    volScalar:Internal    :", os, objects
        );
        print<volVectorField::Internal>
        (
            "    volVector:Internal    :", os, objects
        );
        print<volSphericalTensorField::Internal>
        (
            "    volSphTensor:Internal :", os, objects
        );
        print<volSymmTensorField::Internal>
        (
            "    volSymTensor:Internal :", os, objects
        );
        print<volTensorField::Internal>
        (
            "    volTensor:Internal    :", os, objects
        );
    }


    //- Supported point field types
    static void point(Ostream& os, const IOobjectList& objects)
    {
    }


    //- Supported area field types
    static void area(Ostream& os, const IOobjectList& objects)
    {
        print<areaScalarField>
        (
            "    areaScalar    :", os, objects
        );
        print<areaVectorField>
        (
            "    areaVector    :", os, objects
        );
        print<areaSphericalTensorField>
        (
            "    areaSphTensor :", os, objects
        );
        print<areaSymmTensorField>
        (
            "    areaSymTensor :", os, objects
        );
        print<areaTensorField>
        (
            "    areaTensor    :", os, objects
        );
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

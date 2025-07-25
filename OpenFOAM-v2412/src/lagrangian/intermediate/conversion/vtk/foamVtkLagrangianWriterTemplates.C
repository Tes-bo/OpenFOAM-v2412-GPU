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

\*---------------------------------------------------------------------------*/

#include "foamVtkLagrangianWriter.H"
#include "IOobjectList.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::lagrangianWriter::write(const IOField<Type>& field)
{
    // Other integral types (eg, bool etc) would need cast/convert to label.
    // Similarly for labelVector etc.

    // Ensure consistent output width
    // for (const Type& val : field)
    // {
    //     for (int cmpt=0; cmpt < nCmpt; ++cmpt)
    //     {
    //         format().write(label(component(val, cmpt)));
    //     }
    // }

    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
        vtk::fileWriter::writeBasicField<Type>(field.name(), field);
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
        vtk::fileWriter::writeBasicField<Type>(field.name(), field);
    }
    else
    {
        reportBadState
        (
            FatalErrorInFunction,
            outputState::CELL_DATA,
            outputState::POINT_DATA
        )
            << " for field " << field.name() << nl << endl
            << exit(FatalError);
    }
}


template<class Type>
Foam::label Foam::vtk::lagrangianWriter::writeFields
(
    const wordList& fieldNames,
    bool verbose
)
{
    const fileName localDir(cloudDir());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        // Globally the field is expected to exist (MUST_READ), but can
        // be missing on a local processor.
        //
        // However, constructing IOField with MUST_READ and valid=false fails.
        // Workaround: READ_IF_PRESENT and verify the header globally

        IOobject io
        (
            fieldName,
            mesh_.time().timeName(),
            localDir,
            mesh_,
            IOobject::READ_IF_PRESENT
        );

        // Check global existence to avoid any errors
        if (!returnReduceOr(io.typeHeaderOk<IOField<Type>>(false)))
        {
            continue;
        }

        if (verbose)
        {
            if (!nFields)
            {
                Info<< "        " << pTraits<Type>::typeName << ":";
            }

            Info<< " " << fieldName;
        }

        IOField<Type> field(io);
        this->write<Type>(field);

        ++nFields;
    }

    if (verbose && nFields)
    {
        Info << endl;
    }

    return nFields;
}


template<class Type>
Foam::label Foam::vtk::lagrangianWriter::writeFields
(
    const IOobjectList& objects,
    const bool verbose
)
{
    return writeFields<Type>
    (
        objects.allNames<IOField<Type>>(),  // These are in sorted order
        verbose
    );
}


// ************************************************************************* //

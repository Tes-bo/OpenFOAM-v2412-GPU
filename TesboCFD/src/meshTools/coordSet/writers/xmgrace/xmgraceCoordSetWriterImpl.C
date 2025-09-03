/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "IOmanip.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::coordSetWriters::xmgraceWriter::writeTemplate
(
    const word& fieldName,
    const UPtrList<const Field<Type>>& fieldPtrs
)
{
    if (coords_.size() != fieldPtrs.size())
    {
        FatalErrorInFunction
            << "Attempted to write field: " << fieldName
            << " (" << fieldPtrs.size() << " entries) for "
            << coords_.size() << " sets" << nl
            << exit(FatalError);
    }


    // Field:  rootdir/<TIME>/<field>_setName.agr

    fileName outputFile = getFieldPrefixedPath(fieldName, "agr");

    if (verbose_)
    {
        Info<< "Writing field " << fieldName;
        Info<< " to " << outputFile << endl;
    }

    // Master only
    {
        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os(outputFile, streamOpt_);
        os.precision(precision_);

        // Preamble
        {
            const coordSet& coords = coords_[0];

            os  << "@g0 on" << nl
                << "@with g0" << nl
                << "@    title \"" << coords.name() << '"' << nl
                << "@    xaxis label \"" << coords.axis() << '"' << nl;
        }

        const label setNumber = 0;

        // Plot entry
        os  << "@    s" << setNumber
            << " legend \"" << fieldName << '"' << nl
            << "@target G0.S" << setNumber << nl;

        forAll(coords_, tracki)
        {
            writeTable(os, coords_[tracki], fieldPtrs[tracki], " \t");
        }

        os  << '&' << nl;
        os  << "# end_data" << nl;
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

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

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Output coordinate header
    static void writeCoordHeader
    (
        Ostream& os,
        const coordSet& coords,
        const label count  /* (future?) */
    )
    {
        if (coords.hasVectorAxis())
        {
            os << "x,y,z";
        }
        else
        {
            // axis name
            // const word axisName(coords.axis());
            os << word(coords.axis());
        }
    }

    // Write field name, use named components for VectorSpace
    template<class Type>
    static inline void writeHeader(Ostream& os, const word& fieldName)
    {
        /// os  << ' ';  // Extra space

        const auto nCmpts(pTraits<Type>::nComponents);

        if (pTraits<Type>::rank || nCmpts > 1)
        {
            for (direction d = 0; d < nCmpts; ++d)
            {
                os  << ',' << fieldName
                    << '_' << pTraits<Type>::componentNames[d];
            }
        }
        else
        {
            os  << ',' << fieldName;
        }
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::coordSetWriters::csvWriter::writeTemplate
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

    label nPoints = 0;
    for (const auto& pts : coords_)
    {
        nPoints += pts.size();
    }


    // Field:  rootdir/<TIME>/<field>_setName.csv

    fileName outputFile = getFieldPrefixedPath(fieldName, "csv");

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

        // Header
        {
            writeCoordHeader(os, coords_[0], nPoints);
            writeHeader<Type>(os, fieldName);
            os << nl;
        }

        forAll(coords_, tracki)
        {
            writeTable(os, coords_[tracki], fieldPtrs[tracki], ",");
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "IOmanip.H"
#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Emit x,y,z
    static inline void writePoint(Ostream& os, const point& p)
    {
        os << p.x() << ' ' << p.y() << ' ' << p.z();
    }

    // Emit each component
    template<class Type>
    static inline void writeData(Ostream& os, const Type& val)
    {
        for (direction d = 0; d < pTraits<Type>::nComponents; ++d)
        {
            os  << ' ' << component(val, d);
        }
    }

    // Write area header
    static inline void writeHeaderArea(Ostream& os)
    {
        os  << "  area_x area_y area_z";
    }

    // Write field name, use named components for VectorSpace
    template<class Type>
    static inline void writeHeader(Ostream& os, const word& fieldName)
    {
        os  << ' ';

        const auto nCmpts(pTraits<Type>::nComponents);

        if (pTraits<Type>::rank || nCmpts > 1)
        {
            for (direction d = 0; d < nCmpts; ++d)
            {
                os  << ' ' << fieldName
                    << '_' << pTraits<Type>::componentNames[d];
            }
        }
        else
        {
            os  << ' ' << fieldName;
        }
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::surfaceWriters::rawWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& localValues
)
{
    checkOpen();

    // Field:  rootdir/<TIME>/<field>_surfaceName.raw

    fileName outputFile = outputPath_.path();
    if (useTimeDir() && !timeName().empty())
    {
        // Splice in time-directory
        outputFile /= timeName();
    }

    // Append <field>_surfaceName.raw
    outputFile /= fieldName + '_' + outputPath_.name();
    outputFile.ext("raw");


    // Implicit geometry merge()
    tmp<Field<Type>> tfield = adjustField(fieldName, mergeField(localValues));

    if (verbose_)
    {
        Info<< " to " << outputFile << endl;
    }

    // const meshedSurf& surf = surface();
    const meshedSurfRef& surf = adjustSurface();

    if (UPstream::master() || !parallel_)
    {
        const auto& values = tfield();
        const pointField& points = surf.points();
        const faceList& faces = surf.faces();
        const bool withFaceNormal = (writeNormal_ && !this->isPointData());

        if (!isDir(outputFile.path()))
        {
            mkDir(outputFile.path());
        }

        OFstream os(outputFile, streamOpt_);
        os.precision(precision_);

        // Header
        {
            os  << "# " << fieldName;
            if (this->isPointData())
            {
                os  << "  POINT_DATA ";
            }
            else
            {
                os << "  FACE_DATA ";
            }
            os << values.size() << nl;

            os  << "# x y z";
            writeHeader<Type>(os, fieldName);
            if (withFaceNormal)
            {
                writeHeaderArea(os);
            }
            os  << nl;
        }


        if (this->isPointData())
        {
            // Node values
            forAll(values, elemi)
            {
                writePoint(os, points[elemi]);
                writeData(os, values[elemi]);
                os << nl;
            }
        }
        else
        {
            // Face values
            forAll(values, elemi)
            {
                const face& f = faces[elemi];

                writePoint(os, f.centre(points));
                writeData(os,  values[elemi]);
                if (withFaceNormal)
                {
                    os << ' ';
                    writePoint(os, f.areaNormal(points));
                }
                os << nl;
            }
        }
    }

    wroteGeom_ = true;
    return outputFile;
}


// ************************************************************************* //

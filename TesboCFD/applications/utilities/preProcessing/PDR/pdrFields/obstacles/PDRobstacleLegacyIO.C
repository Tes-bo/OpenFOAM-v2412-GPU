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

\*---------------------------------------------------------------------------*/

#include "PDRobstacle.H"
#include "vector.H"
#include "stringOps.H"
#include "unitConversion.H"
#include <cmath>

#define ReportLineInfo(line, file)                                            \
    if (line >= 0 && !file.empty())                                           \
    {                                                                         \
        Info<< " Line " << line << " of file '" << file << '\'';              \
    }                                                                         \
    Info<< nl;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::PDRobstacle::setFromLegacy
(
    const int groupTypeId,
    const string& buffer,
    const label lineNo,
    const word& inputFile
)
{
    // Handling the optional identifier string can be a pain.
    // Generally can only say that it exists if there are 15 or
    // more columns.
    //
    // Cylinder has 8 normal entries
    // Cuboid, diagonal beam etc have 14 normal entries
    // However, reject anything that looks like a slipped numeric

    double dummy1;

    string in_ident;

    const auto columns = stringOps::splitSpace(buffer);

    for (std::size_t coli = 14; coli < columns.size(); ++coli)
    {
        // See if it can parse into a numerical value
        if (!readDouble(columns[coli].str(), dummy1))
        {
            // Not a numeric value. This must be our identifier
            in_ident = buffer.substr(columns[coli].first - buffer.begin());

            #ifdef FULLDEBUG
            Info<< "Identifier: " << in_ident << nl;
            #endif
            break;
        }
    }

    // Strip off group number
    groupId = groupTypeId / 100;
    typeId  = groupTypeId % 100;

    // This is a safe value
    orient = vector::X;

    switch (typeId)
    {
        case PDRobstacle::CYLINDER:
        {
            // 8 Tokens
            // "%d %lf %lf %lf %lf %lf %d %lf"
            // USP 13/8/14  Read vbkge in case a negative cyl to punch a circular hole

            int in_typeId;
            double in_x, in_y, in_z;
            double in_len, in_dia;
            int in_orient;
            double in_poro;

            int nread =
                sscanf
                (
                    buffer.c_str(),
                    "%d %lf %lf %lf %lf %lf %d %lf",
                    &in_typeId, &in_x, &in_y, &in_z,
                    &in_len, &in_dia, &in_orient,
                    &in_poro
                );

            if (nread < 8)
            {
                Info<< "Expected 8 items, but read in " << nread;
                ReportLineInfo(lineNo, inputFile);
            }

            identifier = in_ident;
            pt = point(in_x, in_y, in_z);

            len() = in_len;
            dia() = in_dia;

            orient = vector::X;  // Check again later

            // Read porosity. Convert to blockage.
            vbkge = 1.0 - in_poro;

            // Orientation (1,2,3) on input -> (0,1,2)
            // - sortBias for later position sorting
            switch (in_orient)
            {
                case 1:
                    orient = vector::X;
                    sortBias = len();
                    break;
                case 2:
                    orient = vector::Y;
                    sortBias = 0.5*dia();
                    break;
                case 3:
                    orient = vector::Z;
                    sortBias = 0.5*dia();
                    break;
                default:
                    sortBias = len();
                    Info<< "Unexpected orientation " << in_orient;
                    ReportLineInfo(lineNo, inputFile);
                    break;
            }
        }
        break;

        case PDRobstacle::DIAG_BEAM:
        {
            // A diagonal block

            // 14 columns + identifier
            // "%d %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %d %lf %s"
            // vbkge (porosity at this stage) should be 0. Not used (yet)

            int in_typeId;
            double in_x, in_y, in_z;
            double in_len, in_theta;
            int in_orient;
            double in_wa, in_wb, in_poro;
            double col_11, col_12, col_14;
            int col_13;

            int nread =
                sscanf
                (
                    buffer.c_str(),
                    "%d %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %d %lf",
                    &in_typeId, &in_x, &in_y, &in_z,
                    &in_len, &in_theta, &in_orient,
                    &in_wa, &in_wb, &in_poro,
                    &col_11, &col_12, &col_13, &col_14
                );

            if (nread < 14)
            {
                Info<< "Expected min 10 items, but read in " << nread;
                ReportLineInfo(lineNo, inputFile);
            }

            identifier = in_ident;
            pt = point(in_x, in_y, in_z);

            len() = in_len;
            dia() = 0;
            theta() = 0; // Fix later on

            orient = vector::X;  // Check again later

            wa  = in_wa;
            wb  = in_wb;

            // Degrees on input, limit to range [0, PI]
            while (in_theta > 180) in_theta -= 180;
            while (in_theta < 0) in_theta += 180;

            // Swap axes when theta > PI/2
            // For 89-90 degrees it becomes -ve, which is picked up
            // in next section
            if (in_theta > 89)
            {
                in_theta -= 90;
                // Swap wa <-> wb
                wa = in_wb;
                wb = in_wa;
            }

            theta() = degToRad(in_theta);

            // Orientation (1,2,3) on input -> (0,1,2)
            // - sortBias for later position sorting
            switch (in_orient)
            {
                case 1:
                    orient = vector::X;
                    sortBias = len();
                    break;

                case 2:
                    orient = vector::Y;
                    sortBias = 0.5*(wa * sin(theta()) + wb * cos(theta()));
                    break;

                case 3:
                    orient = vector::Z;
                    sortBias = 0.5*(wa * cos(theta()) + wb * sin(theta()));
                    break;

                default:
                    sortBias = len();
                    Info<< "Unexpected orientation " << in_orient;
                    ReportLineInfo(lineNo, inputFile);
                    break;
            }


            // If very nearly aligned with axis, turn it into normal block,
            // to avoid 1/tan(theta) blowing up
            if (in_theta < 1)
            {
                switch (orient)
                {
                    case vector::X:
                        span = vector(len(), wa, wb);
                        // Was end center, now lower corner
                        pt.y() = pt.y() - 0.5 * span.y();
                        pt.z() = pt.z() - 0.5 * span.z();
                        break;

                    case vector::Y:
                        span = vector(wb, len(), wa);
                        // Was end center, now lower corner
                        pt.z() = pt.z() - 0.5 * span.z();
                        pt.x() = pt.x() - 0.5 * span.x();
                        break;

                    case vector::Z:
                        span = vector(wa, wb, len());
                        // Was end center, now lower corner
                        pt.x() = pt.x() - 0.5 * span.x();
                        pt.y() = pt.y() - 0.5 * span.y();
                        break;
                }

                typeId = PDRobstacle::CUBOID;
                sortBias = 0;
                xbkge = ybkge = zbkge = vbkge = 1.0;
                blowoff_type = 0;

                Info<< "... changed to type cuboid" << nl;
                break;
            }
        }
        break;

        case PDRobstacle::CUBOID_1:       // Cuboid "Type 1"
        case PDRobstacle::LOUVRE_BLOWOFF: // Louvred wall or blow-off panel
        case PDRobstacle::CUBOID:         // Cuboid
        case PDRobstacle::WALL_BEAM:      // Beam against wall (treated here as normal cuboid)
        case PDRobstacle::GRATING:        // Grating
        case PDRobstacle::RECT_PATCH:     // Inlet, outlet, other b.c. (rectangular)
        {
            // 14 columns + identifier
            // "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf"

            int in_typeId;
            double in_x, in_y, in_z;
            double in_delx, in_dely, in_delz;
            double in_poro, in_porox, in_poroy, in_poroz;
            double col_12;
            int col_13;
            double in_blowoff_time = 0;

            int nread =
                sscanf
                (
                    buffer.c_str(),
                    "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf",
                    &in_typeId, &in_x, &in_y, &in_z,
                    &in_delx, &in_dely, &in_delz,
                    &in_poro, &in_porox, &in_poroy, &in_poroz,
                    &col_12, &col_13, &in_blowoff_time
                );

            blowoff_time = scalar(in_blowoff_time);

            if (nread < 14)
            {
                Info<< "Expected 14 items, but read in " << nread;
                ReportLineInfo(lineNo, inputFile);
            }

            identifier = in_ident;
            pt = point(in_x, in_y, in_z);

            span = vector(in_delx, in_dely, in_delz);

            // Read porosity. Convert to blockage.
            vbkge = 1.0 - in_poro;
            xbkge = 1.0 - in_porox;
            ybkge = 1.0 - in_poroy;
            zbkge = 1.0 - in_poroz;

            if
            (
                typeId == PDRobstacle::CUBOID_1
             || typeId == PDRobstacle::WALL_BEAM
             || typeId == PDRobstacle::RECT_PATCH
            )
            {
                // Check for invalid input

                if (vbkge != 1.0 || xbkge != 1.0 || ybkge != 1.0 || zbkge != 1.0)
                {
                    Info<< "Type " << typeId << " is porous (setting to blockage).";
                    ReportLineInfo(lineNo, inputFile);

                    vbkge = 1;
                    xbkge = 1;
                    ybkge = 1;
                    zbkge = 1;
                }
                if (typeId == PDRobstacle::RECT_PATCH)
                {
                    // Correct interpretation of column 13
                    inlet_dirn = col_13;

                    if (identifier.empty())
                    {
                        FatalErrorInFunction
                            << "RECT_PATCH without a patch name"
                            << exit(FatalError);
                    }
                }
            }
            else if (typeId == PDRobstacle::CUBOID)
            {
            }
            else
            {
                if (!equal(cmptProduct(span), 0))
                {
                    Info<< "Type " << typeId << " has non-zero thickness.";
                    ReportLineInfo(lineNo, inputFile);
                }
            }

            if (typeId == PDRobstacle::LOUVRE_BLOWOFF)
            {
                // Blowoff panel
                blowoff_press = barToPa(col_12);
                blowoff_type = col_13;

                if (blowoff_type == 1)
                {
                    Info<< "Type " << typeId
                        << ": blowoff-type 1 not yet implemented.";
                    ReportLineInfo(lineNo, inputFile);

                    if (blowoff_time != 0)
                    {
                        Info<< "Type " << typeId << " has blowoff time set,"
                            << " not set to blow off cell-by-cell";
                        ReportLineInfo(lineNo, inputFile);
                    }
                }
                else
                {
                    if
                    (
                        (blowoff_type == 1 || blowoff_type == 2)
                     && (col_12 > 0)
                    )
                    {
                        if (col_12 > maxBlowoffPressure)
                        {
                            Info<< "Blowoff pressure (" << col_12
                                << ") too high for blowoff type "
                                << blowoff_type;
                            ReportLineInfo(lineNo, inputFile);
                        }
                    }
                    else
                    {
                        Info<< "Problem with blowoff parameters";
                        ReportLineInfo(lineNo, inputFile);
                        Info<< "Col12 " << col_12
                            << " Blowoff type " << blowoff_type
                            << ", blowoff pressure " << blowoff_press << nl;
                    }
                }
            }
            else if (typeId == PDRobstacle::WALL_BEAM)
            {
                // WALL_BEAM against walls only contribute half to drag
                // if ((col_12 == 1) || (col_12 == -1)) { against_wall_fac = 0.5; }
            }
            else if (typeId == PDRobstacle::GRATING)
            {
                if (col_12 > 0)
                {
                    slat_width = col_12;
                }
                else
                {
                    slat_width = 0;
                }

                // Set orientation
                if (equal(span.x(), 0))
                {
                    orient = vector::X;
                }
                else if (equal(span.y(), 0))
                {
                    orient = vector::Y;
                }
                else
                {
                    orient = vector::Z;
                }
            }
        }
        break;

        case 0:  // Group location
        case PDRobstacle::OLD_INLET: // Ventilation source only
            return false;
            break;

        case PDRobstacle::IGNITION:  // Ignition (now ignored. 2019-04)
            Info<< "Ignition cell type ignored";
            ReportLineInfo(lineNo, inputFile);
            return false;
            break;

        default:
            Info<< "Unexpected type " << typeId;
            ReportLineInfo(lineNo, inputFile);
            return false;
            break;
    }

    return true;
}


// ************************************************************************* //

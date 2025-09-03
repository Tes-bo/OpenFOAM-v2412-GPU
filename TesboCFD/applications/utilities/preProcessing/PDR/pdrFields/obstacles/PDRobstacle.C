/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
#include "boundBox.H"
#include "meshedSurf.H"
#include "axisAngleRotation.H"
#include "coordinateSystem.H"
#include "foamVtkSurfaceWriter.H"
#include "unitConversion.H"
#include "addToMemberFunctionSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineMemberFunctionSelectionTable(PDRobstacle, read, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRobstacle::PDRobstacle()
:
    groupId(0),
    typeId(0),
    orient(vector::X),
    sortBias(0),
    pt(Zero),
    span(Zero),
    wa(0),
    wb(0),
    vbkge(0),
    xbkge(0),
    ybkge(0),
    zbkge(0),
    blowoff_type(0),
    identifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRobstacle::clear()
{
    groupId = 0;
    typeId = 0;
    orient = vector::X;
    sortBias = 0;
    pt = Zero;
    span = Zero;
    wa = 0;
    wb = 0;
    vbkge = 0;
    xbkge = 0;
    ybkge = 0;
    zbkge = 0;
    blowoff_type = 0;
    identifier.clear();
}


void Foam::PDRobstacle::readProperties(const dictionary& dict)
{
    PDRobstacle::clear();

    // Read as word, which handles quoted or unquoted entries
    word obsName;

    if (dict.readIfPresent("name", obsName))
    {
        identifier = std::move(obsName);
    }
}


void Foam::PDRobstacle::scale(const scalar factor)
{
    if (factor <= 0)
    {
        return;
    }

    sortBias *= factor;

    switch (typeId)
    {
        case PDRobstacle::CYLINDER:
        {
            pt *= factor;

            dia() *= factor;
            len() *= factor;
            break;
        }

        case PDRobstacle::DIAG_BEAM:
        {
            pt *= factor;

            len() *= factor;
            wa *= factor;
            wb *= factor;
            break;
        }

        case PDRobstacle::CUBOID_1:
        case PDRobstacle::LOUVRE_BLOWOFF:
        case PDRobstacle::CUBOID:
        case PDRobstacle::WALL_BEAM:
        case PDRobstacle::GRATING:
        case PDRobstacle::RECT_PATCH:
        {
            pt *= factor;
            span *= factor;

            if (typeId == PDRobstacle::GRATING)
            {
                slat_width *= factor;
            }
            break;
        }
    }
}


Foam::scalar Foam::PDRobstacle::volume() const
{
    scalar vol = 0;

    switch (typeId)
    {
        case PDRobstacle::CYLINDER:
            vol = 0.25 * mathematical::pi * sqr(dia()) * len();
            break;

        case PDRobstacle::DIAG_BEAM:
            vol = wa * wb * len();
            break;

        case PDRobstacle::CUBOID_1:
        case PDRobstacle::LOUVRE_BLOWOFF:
        case PDRobstacle::CUBOID:
        case PDRobstacle::WALL_BEAM:
        case PDRobstacle::GRATING:
            vol = cmptProduct(span) * vbkge;
            break;
    }

    return vol;
}


bool Foam::PDRobstacle::tooSmall(const scalar minWidth) const
{
    if (minWidth <= 0)
    {
        return false;
    }

    switch (typeId)
    {
        case PDRobstacle::CYLINDER:
        {
            // The effective half-width
            if ((0.25 * dia() * sqrt(mathematical::pi)) <= minWidth)
            {
                return true;
            }
            break;
        }

        case PDRobstacle::DIAG_BEAM:
        {
            if
            (
                (len() <= minWidth && wa <= minWidth)
             || (len() <= minWidth && wb <= minWidth)
             || (wa  <= minWidth && wb <= minWidth)
            )
            {
                return true;
            }
            break;
        }

        case PDRobstacle::CUBOID_1:
        case PDRobstacle::LOUVRE_BLOWOFF:
        case PDRobstacle::CUBOID:
        case PDRobstacle::WALL_BEAM:
        case PDRobstacle::GRATING:
        case PDRobstacle::RECT_PATCH:
        {
            if
            (
                (span.x() <= minWidth && span.y() <= minWidth)
             || (span.y() <= minWidth && span.z() <= minWidth)
             || (span.z() <= minWidth && span.x() <= minWidth)
            )
            {
                return true;
            }

            break;
        }
    }

    return false;
}


Foam::volumeType Foam::PDRobstacle::trim(const boundBox& bb)
{
    volumeType::type vt = volumeType::UNKNOWN;

    if (!bb.good() || !typeId)
    {
        return vt;
    }

    switch (typeId)
    {
        case PDRobstacle::CYLINDER:
        {
            const scalar rad = 0.5*dia();

            direction e1 = vector::X;
            direction e2 = vector::Y;
            direction e3 = vector::Z;

            if (orient == vector::X)
            {
                e1 = vector::Y;
                e2 = vector::Z;
                e3 = vector::X;
            }
            else if (orient == vector::Y)
            {
                e1 = vector::Z;
                e2 = vector::X;
                e3 = vector::Y;
            }
            else
            {
                orient = vector::Z;  // extra safety?
            }

            if
            (
                (pt[e1] + rad <= bb.min()[e1])
             || (pt[e2] + rad <= bb.min()[e2])
             || (pt[e3] + len() <= bb.min()[e3])
             || (pt[e1] - rad >= bb.max()[e1])
             || (pt[e2] - rad >= bb.max()[e2])
             || (pt[e3] >= bb.max()[e3])
            )
            {
                // No overlap
                return volumeType::OUTSIDE;
            }

            vt = volumeType::INSIDE;

            // Trim obstacle length as required
            if (pt[e3] < bb.min()[e3])
            {
                vt = volumeType::MIXED;
                len() -= bb.min()[e3] - pt[e3];
                pt[e3] = bb.min()[e3];
            }

            if (pt[e3] + len() > bb.max()[e3])
            {
                vt = volumeType::MIXED;
                len() = bb.max()[e3] - pt[e3];
            }

            // Cannot trim diameter very well, so just mark as protruding
            if
            (
                (pt[e1] - rad < bb.min()[e1]) || (pt[e1] + rad > bb.max()[e1])
             || (pt[e2] - rad < bb.min()[e2]) || (pt[e2] + rad > bb.max()[e2])
            )
            {
                vt = volumeType::MIXED;
            }

            break;
        }


        case PDRobstacle::DIAG_BEAM:
        {
            // Not implemented
            break;
        }


        case PDRobstacle::CUBOID_1:
        case PDRobstacle::LOUVRE_BLOWOFF:
        case PDRobstacle::CUBOID:
        case PDRobstacle::WALL_BEAM:
        case PDRobstacle::GRATING:
        case PDRobstacle::RECT_PATCH:
        {
            for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
            {
                if
                (
                    ((pt[cmpt] + span[cmpt]) < bb.min()[cmpt])
                 || (pt[cmpt] > bb.max()[cmpt])
                )
                {
                    // No overlap
                    return volumeType::OUTSIDE;
                }
            }


            vt = volumeType::INSIDE;

            // Trim obstacle as required

            for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
            {
                if (pt[cmpt] < bb.min()[cmpt])
                {
                    vt = volumeType::MIXED;
                    if (span[cmpt] > 0)
                    {
                        span[cmpt] -= bb.min()[cmpt] - pt[cmpt];
                    }
                    pt[cmpt] = bb.min()[cmpt];
                }


                if (pt[cmpt] + span[cmpt] > bb.max()[cmpt])
                {
                    vt = volumeType::MIXED;
                    span[cmpt] -= bb.max()[cmpt] - pt[cmpt];
                }
            }

            break;
        }
    }

    return vt;
}


Foam::meshedSurface Foam::PDRobstacle::surface() const
{
    meshedSurface surf;

    const PDRobstacle& obs = *this;

    switch (obs.typeId)
    {
        case PDRobstacle::CUBOID_1 :
        case PDRobstacle::CUBOID :
        {
            boundBox box(obs.pt, obs.pt + obs.span);

            pointField pts(box.points());
            faceList fcs(boundBox::hexFaces());

            surf.transfer(pts, fcs);

            break;
        }

        case PDRobstacle::DIAG_BEAM :
        {
            boundBox box(Zero);

            switch (orient)
            {
                case vector::X:
                {
                    box.min() = vector(0, -0.5*obs.wa, -0.5*obs.wb);
                    box.max() = vector(obs.len(), 0.5*obs.wa,  0.5*obs.wb);
                    break;
                }

                case vector::Y:
                {
                    box.min() = vector(-0.5*obs.wb, 0, -0.5*obs.wa);
                    box.max() = vector(0.5*obs.wb, obs.len(), 0.5*obs.wa);
                    break;
                }

                case vector::Z:
                {
                    box.min() = vector(-0.5*obs.wa, -0.5*obs.wb, 0);
                    box.max() = vector(0.5*obs.wa,  0.5*obs.wb, obs.len());
                    break;
                }
            }

            coordinateSystem cs
            (
                obs.pt,
                coordinateRotations::axisAngle
                (
                    vector::components(obs.orient),
                    obs.theta(),
                    false
                )
            );

            pointField pts0(box.points());
            faceList fcs(boundBox::hexFaces());

            pointField pts(cs.globalPosition(pts0));

            surf.transfer(pts, fcs);

            break;
        }

        case PDRobstacle::CYLINDER :
        {
            // Tessellation 12 looks fairly reasonable

            constexpr int nDiv = 12;

            point org(obs.pt);

            direction e1 = vector::X;
            direction e2 = vector::Y;
            direction e3 = vector::Z;

            if (obs.orient == vector::X)
            {
                e1 = vector::Y;
                e2 = vector::Z;
                e3 = vector::X;
            }
            else if (obs.orient == vector::Y)
            {
                e1 = vector::Z;
                e2 = vector::X;
                e3 = vector::Y;
            }

            pointField pts(2*nDiv, org);
            faceList   fcs(2 + nDiv);

            // Origin for back
            org[e3] += obs.len();
            SubList<point>(pts, nDiv, nDiv) = org;

            const scalar radius = 0.5*obs.dia();

            for (label i=0; i < nDiv; ++i)
            {
                const scalar angle = (i * mathematical::twoPi) / nDiv;
                const scalar s = radius * sin(angle);
                const scalar c = radius * cos(angle);

                pts[i][e1] += s;
                pts[i][e2] += c;

                pts[nDiv+i][e1] += s;
                pts[nDiv+i][e2] += c;
            }

            // Side-faces
            for (label facei=0; facei < nDiv; ++facei)
            {
                face& f = fcs[facei];
                f.resize(4);

                f[0] = facei;
                f[3] = (facei + 1) % nDiv;
                f[1] = f[0] + nDiv;
                f[2] = f[3] + nDiv;
            }

            {
                // Front face
                face& f1 = fcs[nDiv];
                f1.resize(nDiv);

                f1[0] = 0;
                for (label pti=1; pti < nDiv; ++pti)
                {
                    f1[pti] = nDiv-pti;
                }

                // Back face
                labelList& f2 = fcs[nDiv+1];
                f2 = identity(nDiv, nDiv);
            }

            surf.transfer(pts, fcs);

            break;
        }

        case PDRobstacle::RECT_PATCH :
        {
            pointField pts(4, obs.span);
            pts[0] = Zero;

            switch (obs.inlet_dirn)
            {
                case -1:
                case 1:
                {
                    for (auto& p : pts)
                    {
                        p.x() = 0;
                    }

                    pts[1].z() = 0;
                    pts[3].y() = 0;
                    break;
                }
                case -2:
                case 2:
                {
                    for (auto& p : pts)
                    {
                        p.y() = 0;
                    }

                    pts[1].x() = 0;
                    pts[3].z() = 0;
                    break;
                }
                default:
                {
                    for (auto& p : pts)
                    {
                        p.z() = 0;
                    }

                    pts[1].y() = 0;
                    pts[3].x() = 0;
                    break;
                }
            }

            // pts += obs.pt;

            faceList fcs(one{}, face(identity(4)));

            surf.transfer(pts, fcs);

            break;
        }

        default:
            break;

//         LOUVRE_BLOWOFF = 5,
//         WALL_BEAM =  7,
//         GRATING   =  8,
//         CIRC_PATCH  = 12,
//         MESH_PLANE  = 46,
    }

    return surf;
}


Foam::label Foam::PDRobstacle::addPieces
(
    vtk::surfaceWriter& surfWriter,
    const UList<PDRobstacle>& list,
    label pieceId
)
{
    for (const PDRobstacle& obs : list)
    {
        meshedSurface surf(obs.surface());

        if (!surf.empty())
        {
            surfWriter.piece(surf.points(), surf.surfFaces());

            surfWriter.writeGeometry();
            surfWriter.beginCellData(2);
            surfWriter.writeUniform("group", label(obs.groupId));
            surfWriter.writeUniform("type", label(obs.typeId));
            surfWriter.writeUniform("obstacle", pieceId);
            ++pieceId;
        }
    }

    return pieceId;
}


void Foam::PDRobstacle::generateVtk
(
    const fileName& outputDir,
    const UList<PDRobstacle>& obslist,
    const UList<PDRobstacle>& cyllist
)
{
    label pieceId = 0;

    vtk::surfaceWriter surfWriter
    (
        pointField::null(),
        faceList::null(),
        // vtk::formatType::INLINE_ASCII,
        (outputDir / "Obstacles"),
        false  // serial only
    );

    pieceId = addPieces(surfWriter, obslist, pieceId);
    pieceId = addPieces(surfWriter, cyllist, pieceId);

    Info<< "Wrote " << pieceId << " obstacles (VTK) to "
        << outputDir/"Obstacles" << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<PDRobstacle>& iproxy
)
{
    const auto& obs = *iproxy;

    switch (obs.typeId)
    {
        case PDRobstacle::CUBOID_1 :
        case PDRobstacle::CUBOID :
            os  << "box  { point " << obs.pt
                << "; size " << obs.span
                << "; }";
            break;

        case PDRobstacle::CYLINDER :
            os  << "cyl { point " << obs.pt
                << "; length " << obs.len() << "; diameter " << obs.dia()
                << "; direction " << vector::componentNames[obs.orient]
                << "; }";
            break;

        case PDRobstacle::DIAG_BEAM :
            os  << "diag { point " << obs.pt
                << "; length " << obs.len()
                << "; width (" << obs.wa << ' ' << obs.wb << ')'
                << "; angle " << radToDeg(obs.theta())
                << "; direction " << vector::componentNames[obs.orient]
                << "; }";
            break;

        case PDRobstacle::WALL_BEAM :
            os  << "wallbeam { point " << obs.pt
                << " size " << obs.span
                << "; }";
            break;

        case PDRobstacle::GRATING :
            os  << "grate { point " << obs.pt
                << "; size " << obs.span
                << "; slats " << obs.slat_width
                << "; }";
            break;

        case PDRobstacle::LOUVER_BLOWOFF :
            os  << "louver { point " << obs.pt
                << "; size " << obs.span
                << "; pressure " << paToBar(obs.blowoff_press)
                << "; }";
            break;

        case PDRobstacle::RECT_PATCH :
            os  << "patch { " << obs.pt
                << "; size " << obs.span
                << "; name " << obs.identifier
                << "; }";
            break;

        case PDRobstacle::OLD_INLET :
        case PDRobstacle::OLD_BLOWOFF :
        case PDRobstacle::IGNITION :
            os  << "/* ignored: " << obs.typeId << " */";
            break;

        default:
            os  << "/* unknown: " << obs.typeId << " */";
            break;

    }

    return os;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

bool Foam::operator<(const PDRobstacle& a, const PDRobstacle& b)
{
    return (a.pt.x() + a.sortBias) < (b.pt.x() + b.sortBias);
}


// ************************************************************************* //

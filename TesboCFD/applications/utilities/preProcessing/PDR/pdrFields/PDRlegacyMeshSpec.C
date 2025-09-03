/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
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

// Read spec that resemble this:
//
//   xmesh
//   (
//       ( -8.18  14  0.83 )
//       ( -0.69  11  1.00 )
//       (  0.69  14  1.20 )
//       (  8.18  0 )
//   )
//
//   ymesh
//   (
//       ...
//   )

#include "PDRlegacy.H"

// OpenFOAM includes
#include "error.H"
#include "IFstream.H"
#include "stringOps.H"
#include "OSspecific.H"

#define XMESH_TAG "xmesh"
#define YMESH_TAG "ymesh"
#define ZMESH_TAG "zmesh"


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace PDRlegacy
{
namespace Detail
{

//- Simple structure for processing a mesh spec entry
//  Handles an entry with (float), (float int), or (float int float)
//
//  This is for handling a sub-spec that resembles this:
//
//  (
//      ( -8.18  14  0.83 )
//      ( -0.69  11  1.00 )
//      (  0.69  14  1.20 )
//      (  8.18  0 )
//  )
struct pdrMeshSpecLine
{
    scalar knot;
    label  ndiv;
    scalar factor;

    pdrMeshSpecLine() : knot(0), ndiv(0), factor(0) {}

    //- Cheap means to avoid uniform list output
    bool operator!=(const pdrMeshSpecLine&) const
    {
        return false;
    }
};


// Read mesh-spec entry
Istream& operator>>(Istream& is, pdrMeshSpecLine& spec)
{
    spec.knot = 0;
    spec.ndiv = 0;
    spec.factor = 0;

    is.readBegin("pdrMeshSpecLine");

    // Must have a point
    is >> spec.knot;

    token tok(is);
    if (tok.isLabel())
    {
        spec.ndiv = tok.labelToken();

        if (spec.ndiv)
        {
            is >> spec.factor;
        }
    }
    else
    {
        is.putBack(tok);
    }

    is.readEnd("pdrMeshSpecLine");

    is.check(FUNCTION_NAME);
    return is;
}

#ifdef FULLDEBUG
// Write mesh-spec entry
Ostream& operator<<(Ostream& os, const pdrMeshSpecLine& spec)
{
    os << token::BEGIN_LIST << spec.knot;

    if (spec.ndiv)
    {
        os  << token::SPACE << spec.ndiv
            << token::SPACE << spec.factor;
    }

    os << token::END_LIST;

    return os;
}
#endif


void read_spec(ISstream& is, const direction cmpt, List<scalar>& gridPoint)
{
    if (!gridPoint.empty())
    {
        FatalErrorInFunction
            << "Duplicate specification of "
            << vector::componentNames[cmpt]
            << " grid"
            << exit(FatalError);
    }

    List<pdrMeshSpecLine> specs(is);

    if (specs.size() < 2)
    {
        FatalErrorInFunction
            << "Grid specification for " << vector::componentNames[cmpt]
            << " is too small. Need at least two points!" << nl
            << exit(FatalError);
    }

    specs.last().ndiv = 0; // safety


    DynamicList<scalar> knots;
    DynamicList<label>  divisions;
    DynamicList<scalar> factors;

    for (const auto& spec : specs)
    {
        knots.append(spec.knot);

        if (spec.ndiv < 1)
        {
            break;
        }
        divisions.append(spec.ndiv);
        factors.append(spec.factor);
    }

    label nPoints = 1;
    for (const label nDiv : divisions)
    {
        nPoints += nDiv;
    }

    if (nPoints < 2)
    {
        FatalErrorInFunction
            << "No cells defined for direction "
            << vector::componentNames[cmpt] << nl
            << exit(FatalError);
    }


    // Define the grid points
    gridPoint.resize(nPoints);

    const label nSegments = divisions.size();

    label start = 0;

    for (label segmenti=0; segmenti < nSegments; ++segmenti)
    {
        const label nDiv = divisions[segmenti];
        const scalar factor = factors[segmenti];

        SubList<scalar> subPoint(gridPoint, nDiv+1, start);
        start += nDiv;

        subPoint[0] = knots[segmenti];
        subPoint[nDiv] = knots[segmenti+1];

        const scalar dist = (subPoint.last() - subPoint.first());

        if (equal(factor, scalar(1)))
        {
            for (label i=1; i < nDiv; ++i)
            {
                subPoint[i] = (subPoint[0] + (dist * i)/nDiv);
            }
        }
        else
        {
            scalar delta = dist * (1.0 - factor) / (1.0 - ::pow(factor, nDiv));

            scalar xyz = subPoint[0];

            for (label i=0; i < nDiv; ++i)
            {
                subPoint[i] = xyz;
                xyz += delta;
                delta *= factor;
            }
        }
    }
}


} // End namespace Detail
} // End namespace PDRlegacy
} // End namespace Foam


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::PDRlegacy::print_info(const PDRblock& block)
{
    Info<< "PDRblock" << nl
        << "    nCells: " << block.sizes() << nl
        << "    Box: " << block.bounds() << nl
        << "x " << flatOutput(block.grid().x()) << nl
        << "y " << flatOutput(block.grid().y()) << nl
        << "z " << flatOutput(block.grid().z()) << nl
        << endl;
}


void Foam::PDRlegacy::read_mesh_spec(const fileName& casepath, PDRblock& block)
{
    Info<< "Reading pdrMeshSpec (legacy format)" << nl;

    bool processed = false;

    for (const fileName dirName : { "system", "constant/polyMesh" })
    {
        fileName path
        (
            casepath / dirName / "pdrMeshSpec"
        );

        if (Foam::isFile(path))
        {
            IFstream is(path);

            read_mesh_spec(is, block);
            processed = true;
            break;
        }
    }

    if (!processed)
    {
        FatalErrorInFunction
            << "Did not process pdrMeshSpec" << nl
            << exit(FatalError);
    }
}


void Foam::PDRlegacy::read_mesh_spec(ISstream& is, PDRblock& block)
{
    Vector<scalarList> grid;

    string line;

    while (is.good())
    {
        is.getLine(line);
        stringOps::inplaceTrim(line);

        if (line == XMESH_TAG)
        {
            Detail::read_spec(is, vector::X, grid.x());
        }
        else if (line == YMESH_TAG)
        {
            Detail::read_spec(is, vector::Y, grid.y());
        }
        else if (line == ZMESH_TAG)
        {
            Detail::read_spec(is, vector::Z, grid.z());
        }
    }

    for (direction cmpt=0; cmpt < vector::nComponents; ++cmpt)
    {
        if (grid[cmpt].empty())
        {
            FatalErrorInFunction
                << "No specification for "
                << vector::componentNames[cmpt]
                << " grid" << nl
                << exit(FatalError);
        }
    }

    block.reset(grid.x(), grid.y(), grid.z());

    #ifdef FULLDEBUG
    print_info(block);
    #endif
}


// ************************************************************************* //

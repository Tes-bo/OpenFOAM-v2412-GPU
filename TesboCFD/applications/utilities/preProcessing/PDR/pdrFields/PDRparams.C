/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "PDRparams.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * * Global Data * * * * * * * * * * * * * * * //

// Global parameter settings
Foam::PDRparams Foam::pars;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PDRparams::readDefaults(const dictionary& dict)
{
    dict.readIfPresent("legacyMeshSpec", legacyMeshSpec);
    dict.readIfPresent("legacyObsSpec",  legacyObsSpec);

    dict.readIfPresent("two_d", two_d);
    dict.readIfPresent("yCyclic", yCyclic);
    dict.readIfPresent("ySymmetry", ySymmetry);
    dict.readIfPresent("deluge", deluge);

    dict.readIfPresent("newFields", new_fields);
    dict.readIfPresent("noIntersectN", noIntersectN);

    dict.readIfPresent("blockedFacesWallFn", blockedFacesWallFn);
    dict.readIfPresent("ignoreGratings", ignoreGratings);

    outer_orthog = dict.found("outer_orthog");

    dict.readIfPresent("debug.level", debugLevel);
    dict.readIfPresent("nFacesToBlockC", nFacesToBlockC);
    dict.readIfPresent("nPairsToBlockC", nPairsToBlockC);
    dict.readIfPresent("overlaps", overlaps);

    dict.readIfPresent("gridPointTol", gridPointTol);

    dict.readIfPresent("Cb_r", cb_r);
    dict.readIfPresent("Cb_s", cb_s);

    dict.readIfPresent("Cd_r", cd_r);
    dict.readIfPresent("Cd_s", cd_s);

    dict.readIfPresent("congRegionMaxBetav", cong_max_betav);

    dict.readIfPresent("min_overlap_vol", min_overlap_vol);
    dict.readIfPresent("min_overlap_area", min_overlap_area);
    dict.readIfPresent("min_width", min_width);
    dict.readIfPresent("empty_lobs_fac", empty_lobs_fac);
    dict.readIfPresent("outerCombFac", outerCombFac);
    dict.readIfPresent("obs_expand", obs_expand);

    dict.readIfPresent("def_grating_slat_w", def_grating_slat_w);
    dict.readIfPresent("blockedCellPoros", blockedCellPoros);
    dict.readIfPresent("blockedFacePar", blockedFacePar);
    dict.readIfPresent("maxCR", maxCR);

    dict.readIfPresent("blockageNoCT", blockageNoCT);
    dict.readIfPresent("scale", scale);


    const dictionary* dictptr;

    groundPatchName = "ground";
    outerPatchName = "outer";

    if ((dictptr = dict.findDict("patchNames")) != nullptr)
    {
        const dictionary& d = *dictptr;

        d.readIfPresent("ground", groundPatchName);
        d.readIfPresent("outer", outerPatchName);
    }

    UPatchBc = "fixedValue;value uniform (0 0 0)";
    if (dict.readIfPresent("UPatchBc", UPatchBc))
    {
        stringOps::inplaceTrim(UPatchBc);
    }
}


void Foam::PDRparams::read(const dictionary& dict)
{
    readDefaults(dict);

    dict.readEntry("obsFileDir", obsfile_dir);
    dict.readEntry("obsFileNames", obsfile_names);

    stringOps::inplaceExpand(obsfile_dir);

    for (auto& f : obsfile_names)
    {
        stringOps::inplaceExpand(f);
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::PDRparams

Description
    Parameters for PDRsetFields

SourceFiles
    PDRparams.C

\*---------------------------------------------------------------------------*/

#ifndef PDRparams_H
#define PDRparams_H

#include "labelList.H"
#include "scalarList.H"
#include "wordList.H"
#include "fileNameList.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class PDRparams Declaration
\*---------------------------------------------------------------------------*/

class PDRparams
{
public:

    // Data Members

        fileName obsfile_dir;
        wordList obsfile_names;
        word    timeName;

        //- The name for the "ground" patch
        word    groundPatchName;

        //- The name for the "outer" patch
        word    outerPatchName;

        string  UPatchBc;   //!< "fixedValue;value uniform (0 0 0)"

        bool legacyMeshSpec{false};
        bool legacyObsSpec{false};

        bool two_d{false};
        bool yCyclic{false};
        bool ySymmetry{false};
        bool deluge{false};
        bool new_fields{true};
        bool noIntersectN{true};
        bool blockedFacesWallFn{false};
        bool ignoreGratings{false};
        bool outer_orthog{false};

        int debugLevel{0};

        //- Min number of blocked cell faces
        //- for a cell to be marked as blocked
        int nFacesToBlockC{6};

        //- Min number of blocked cell face pairs (on opposite faces of a cell)
        //- for a cell to be marked as blocked
        int nPairsToBlockC{3};

        //- Flag to control which overlap calculations are performed
        int overlaps{0x7};

        scalar gridPointTol{0.02};

        scalar cb_r{0.035};
        scalar cb_s{0.08};

        scalar cd_r{1.2};
        scalar cd_s{2.0};

        scalar cong_max_betav{1.0};

        scalar min_overlap_vol{0};
        scalar min_overlap_area{0};

        //- Ignore obstacles with second dimension (or diameter) less than this
        scalar min_width{0.001};

        //- Lobs in empty cell is this * cube root of cell volume
        scalar empty_lobs_fac{1.0};

        //- Value for outer region
        scalar outerCombFac{1.0};

        scalar obs_expand{0};

        //- Default slat thickness grating
        scalar def_grating_slat_w{0.005};

        //- Cells with porosity less than this are blocked
        scalar blockedCellPoros{0.05};

        //- Faces with area blockage greater than this are blocked
        scalar blockedFacePar{0.95};

        //- Upper limit on CR (CT also gets limited)
        scalar maxCR{1e30};

        //- If a single obstacle blocks a cell by more than this,
        //- then no CT in that direction
        scalar blockageNoCT{0.95};

        //- Overall scale factor
        scalar scale{1.0};


    // Constructors

        //- Default construct
        PDRparams() = default;


    // Member Functions

        //- Set or read defaults from dictionary.
        //  Can also be used with an empty dictionary
        void readDefaults(const dictionary& dict);

        //- Read program parameters from dictionary
        void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern Foam::PDRparams pars;  //!< Globals for program parameters (ugly hack)

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

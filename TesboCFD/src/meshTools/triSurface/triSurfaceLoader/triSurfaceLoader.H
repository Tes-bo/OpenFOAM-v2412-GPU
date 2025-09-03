/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::triSurfaceLoader

Description
    Convenience class for loading single or multiple surface files
    from the constant/triSurface (or other) directory.

    Surfaces selection based on word, wordRe, list of wordRe.
    If multiple surfaces are selected, they are concatenated into a
    single surface with offset faces,points,regions.

SourceFiles
    triSurfaceLoader.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_triSurfaceLoader_H
#define Foam_triSurfaceLoader_H

#include "triSurface.H"
#include "wordRes.H"
#include "Enum.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class triSurfaceLoader;
class Time;

/*---------------------------------------------------------------------------*\
                      Class triSurfaceLoader Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceLoader
{
public:

    //- The file loading options for triSurfaceLoader
    enum loadingOption
    {
        SINGLE_REGION,  //!< "single" = One region for all files
        FILE_REGION,    //!< "file"   = One region for each file
        OFFSET_REGION,  //!< "offset" = Offset regions per file
        MERGE_REGION    //!< "merge"  = Merge regions by name
    };

    //- The loading enumeration names
    static const Enum<loadingOption> loadingOptionNames;


private:

    // Private Data

        //- The directory to load from (eg, case/constant/triSurface)
        fileName directory_;

        //- All available files
        wordList available_;

        //- Selected files
        wordList selected_;


    // Private Member Functions

        //- No copy construct
        triSurfaceLoader(const triSurfaceLoader&) = delete;

        //- No copy assignment
        void operator=(const triSurfaceLoader&) = delete;

public:

    // Constructors

        //- Construct with directory name
        triSurfaceLoader(const fileName& directory);

        //- Construct with time. Selects "constant/triSurface" directory.
        triSurfaceLoader(const Time& runTime);


    //- Destructor
    ~triSurfaceLoader() = default;


    // Member Functions

      // Access

        //- The directory being used
        const fileName& directory() const noexcept
        {
            return directory_;
        }

        //- The list of available files
        const wordList& available() const noexcept
        {
            return available_;
        }

        //- The list of selected files
        const wordList& selected() const noexcept
        {
            return selected_;
        }


      // Edit

        //- Read directory and populate the 'available' files.
        //  Automatically called on construction.
        label readDir();

        //- Populates 'selected' with all available files.
        label selectAll();

        //- Populates 'selected' with a subset of the available files.
        label select(const word& name);

        //- Populates 'selected' with a subset of the available files.
        label select(const wordRe& mat);

        //- Populates 'selected' with a subset of the available files.
        label select(const UList<wordRe>& matcher);

        //- Load a single file, or load and combine multiple selected files
        //  Optionally scale the surface(s) on input, with a zero or negative
        //  scale factor treated as no scaling.
        autoPtr<triSurface> load
        (
            const enum loadingOption opt = loadingOption::OFFSET_REGION,
            const scalar scaleFactor = -1
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

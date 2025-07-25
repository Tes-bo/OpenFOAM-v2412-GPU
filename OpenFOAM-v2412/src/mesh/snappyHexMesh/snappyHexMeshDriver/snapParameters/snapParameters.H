/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::snapParameters

Description
    Simple container to keep together snap specific information.

SourceFiles
    snapParameters.C

\*---------------------------------------------------------------------------*/

#ifndef snapParameters_H
#define snapParameters_H

#include "dictionary.H"
#include "scalar.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class snapParameters Declaration
\*---------------------------------------------------------------------------*/

class snapParameters
{
    // Private data

        const label nSmoothPatch_;

        const label nSmoothInternal_;

        const scalar snapTol_;

        const label nSmoothDispl_;

        const label nSnap_;

        const label nFeatureSnap_;

        const Switch explicitFeatureSnap_;

        const Switch implicitFeatureSnap_;

        const Switch multiRegionFeatureSnap_;

        const Switch detectNearSurfacesSnap_;

        const Switch strictRegionSnap_;

        const Switch detectBaffles_;

        const Switch baffleFeaturePoints_;

        const Switch releasePoints_;

        const Switch stringFeatures_;

        const Switch avoidDiagonal_;


        //- How often needs face splitting be run
        const label nFaceSplitInterval_;

        //- When is angle too concave too split
        const scalar concaveAngle_;

        //- When is face-split not sufficiently diagonal
        const scalar minAreaRatio_;

        //- in dry-run mode
        const bool dryRun_;


    // Private Member Functions

        //- No copy construct
        snapParameters(const snapParameters&) = delete;

        //- No copy assignment
        void operator=(const snapParameters&) = delete;


public:

    // Constructors

        //- Construct from dictionary
        snapParameters(const dictionary& dict, const bool dryRun = false);


    // Member Functions

        // Access

            //- Number of patch smoothing iterations before finding
            //  correspondence to surface
            label nSmoothPatch() const
            {
                return nSmoothPatch_;
            }

            //- Number of internal point smoothing iterations (combined with
            //  nSmoothPatch
            label nSmoothInternal() const
            {
                return nSmoothInternal_;
            }

            //- Relative distance for points to be attracted by surface
            //  feature point
            //  or edge. True distance is this factor times local
            //  maximum edge length.
            scalar snapTol() const
            {
                return snapTol_;
            }

            //- Number of mesh displacement smoothing iterations.
            label nSmoothDispl() const
            {
                return nSmoothDispl_;
            }

            //- Maximum number of snapping relaxation iterations. Should stop
            //  before upon reaching a correct mesh.
            label nSnap() const
            {
                return nSnap_;
            }


            // Surface snapping specific

                //- Override attraction to nearest with intersection location
                //  at near surfaces
                Switch detectNearSurfacesSnap() const
                {
                    return detectNearSurfacesSnap_;
                }

                //- Attract point to corresponding surface region only
                Switch strictRegionSnap() const
                {
                    return strictRegionSnap_;
                }


            // Feature edge snapping specific

                label nFeatureSnap() const
                {
                    return nFeatureSnap_;
                }

                Switch explicitFeatureSnap() const
                {
                    return explicitFeatureSnap_;
                }

                Switch implicitFeatureSnap() const
                {
                    return implicitFeatureSnap_;
                }

                Switch multiRegionFeatureSnap() const
                {
                    return multiRegionFeatureSnap_;
                }

                Switch detectBaffles() const
                {
                    return detectBaffles_;
                }

                Switch baffleFeaturePoints() const
                {
                    return baffleFeaturePoints_;
                }

                Switch releasePoints() const
                {
                    return releasePoints_;
                }

                Switch stringFeatures() const
                {
                    return stringFeatures_;
                }

                Switch avoidDiagonal() const
                {
                    return avoidDiagonal_;
                }


                // Face splitting

                    label nFaceSplitInterval() const
                    {
                        return nFaceSplitInterval_;
                    }

                    scalar concaveAngle() const
                    {
                        return concaveAngle_;
                    }

                    scalar minAreaRatio() const
                    {
                        return minAreaRatio_;
                    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

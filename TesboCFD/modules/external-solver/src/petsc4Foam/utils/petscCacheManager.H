/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
    Copyright (C) 2020 Simone Bna
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
    Foam::petscCacheManager

Description
    The manager class for the caching of matrix and preconditioner.

    An Example is:

    \verbatim
    petsc
    {
        ...
        caching
        {
            matrix
            {
                update  always;  // (default: always)
            }

            preconditioner
            {
                update  periodic;

                periodicCoeffs
                {
                    frequency   3;
                }
            }
        }
    }
    \endverbatim

SourceFiles
    petscCacheManager.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscUtils_Caching_H
#define Foam_module_petscUtils_Caching_H

#include "Enum.H"
#include "clockValue.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PetscUtils
{

/*---------------------------------------------------------------------------*\
                           Class Caching Declaration
\*---------------------------------------------------------------------------*/

struct Caching
{
    // Public Data Types

        //- Caching update types
        enum updateTypes
        {
            Always,         //!< "always" update every time-step [default]
            Periodic,       //!< "periodic" update at given period
            Adaptive,       //!< "adaptive" update scheme
            Never,          //!< "never" update (or "none")
        };


        //- Names for the update types
        static const Enum<updateTypes> updateTypeNames_;


    // Member Data

        //- Caching update type
        updateTypes updateType_;

        //- Cache update frequency (for Periodic)
        int updateFreq_;

        //- Relative iteration (for Adaptive or Periodic)
        int relIter_;

        //- Elapsed time (s) for current iteration
        double timeIter_;

        //- Elapsed time (s) for first iteration
        double time0_;

        //- Elapsed time (s) for second iteration
        double time1_;

        //- Timer value
        clockValue timer_;


    // Constructors

        //- Default construct. Always update (no caching).
        Caching()
        :
           updateType_{updateTypes::Always},
           updateFreq_{1},
           relIter_{0},
           timeIter_{0},
           time0_{0},
           time1_{0}
        {}


    // Member Functions

        //- The caching type
        updateTypes type() const noexcept
        {
            return updateType_;
        }

        //- True if caching type is Adaptive
        bool isAdaptive() const noexcept
        {
            return updateTypes::Adaptive == updateType_;
        }

        //- True if caching type is Periodic
        bool isPeriodic() const noexcept
        {
            return updateTypes::Periodic == updateType_;
        }


        //- Initialize from dictionay
        void init(const word& key, const dictionary& dict)
        {
            relIter_ = 0;

            dictionary cacheDict = dict.subOrEmptyDict(key);
            updateType_ =
                updateTypeNames_.getOrDefault
                (
                    "update",
                    cacheDict,
                    updateTypes::Always,
                    false  // non-failsafe - Fatal on bad enumeration
                );

            if (updateTypes::Periodic == updateType_)
            {
                dictionary coeffs(cacheDict.subOrEmptyDict("periodicCoeffs"));
                updateFreq_ = coeffs.getOrDefault<int>("frequency", 1);
            }
        }


        //- Begin/restart timing for iterations
        void beginTimes()
        {
            timer_.update();
        }


        //- Update timing for iterations
        void updateTimes()
        {
            // Elapsed timing interval (s)
            timeIter_ = timer_.elapsed();

            if (relIter_ == 0)
            {
                time0_ = timeIter_;
            }
            else if (relIter_ == 1)
            {
                time1_ = timeIter_;
            }
        }


        //- Check if an update is required
        bool needsUpdate()
        {
            // Default: Always update
            bool need = true;

            switch (updateType_)
            {
                case updateTypes::Never:
                {
                    need = false;
                    break;
                }

                case updateTypes::Always:
                {
                    break;
                }

                case updateTypes::Periodic:
                {
                    if
                    (
                        updateFreq_ <= 1
                     || (relIter_ % updateFreq_) == 0
                    )
                    {
                        relIter_ = 0;
                    }
                    else
                    {
                        need = false;
                    }
                    break;
                }

                case updateTypes::Adaptive:
                {
                    if (relIter_ > 3)  // Need at least three iterations
                    {
                        const double ratio0 = time0_ / timeIter_;

                        const double ratio1 = time1_ / timeIter_;

                        const int nsteps =
                            min(1e5, ratio0 * (1. / mag(1. - ratio1 + 1e-6)));

                        if (relIter_ >= nsteps)
                        {
                            relIter_ = 0;
                        }
                        else
                        {
                            need = false;
                        }
                    }
                    break;
                }
            }

            return need;
        }


    // Member Operators

        //- Increment the relative iteration
        void operator++()
        {
            ++relIter_;
        }
};

} // End namespace PetscUtils


/*---------------------------------------------------------------------------*\
                     Class petscCacheManager Declaration
\*---------------------------------------------------------------------------*/

class petscCacheManager
{
    // Private Data

        PetscUtils::Caching matrixCaching;
        PetscUtils::Caching precondCaching;


public:

    // Constructors

        //- Default construct
        petscCacheManager() = default;


    //- Destructor
    ~petscCacheManager() = default;


    void init(const dictionary& dict)
    {
        matrixCaching.init("matrix", dict);
        precondCaching.init("preconditioner", dict);
    }


    // Member Functions

        bool needsMatrixUpdate()
        {
            return matrixCaching.needsUpdate();
        }

        bool needsPrecondUpdate()
        {
            return precondCaching.needsUpdate();
        }

        void eventBegin()
        {
            if (precondCaching.isAdaptive())
            {
                precondCaching.beginTimes();
            }

            // Should we handle/disable Adaptive for matrixCaching?
        }

        void eventEnd()
        {
            if (precondCaching.isAdaptive())
            {
                precondCaching.updateTimes();
                ++precondCaching;
            }
            else if (precondCaching.isPeriodic())
            {
                ++precondCaching;
            }

            // Should we handle/disable Adaptive for matrixCaching?
            if (matrixCaching.isPeriodic())
            {
                ++matrixCaching;
            }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

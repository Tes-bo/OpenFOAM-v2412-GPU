/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::nearestFaceAMI

Description
    Nearest-face Arbitrary Mesh Interface (AMI) method

SourceFiles
    nearestFaceAMI.C

\*---------------------------------------------------------------------------*/

#ifndef nearestFaceAMI_H
#define nearestFaceAMI_H

#include "AMIInterpolation.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class nearestFaceAMI Declaration
\*---------------------------------------------------------------------------*/

class nearestFaceAMI
:
    public AMIInterpolation
{
public:

    typedef Tuple2<pointIndexHit, scalar> nearestAndDist;

    //- Helper class for finding nearest
    class nearestEqOp
    {

    public:

        void operator()(nearestAndDist& x, const nearestAndDist& y) const
        {
            if (y.first().hit())
            {
                if (!x.first().hit())
                {
                    x = y;
                }
                else if (y.second() < x.second())
                {
                    x = y;
                }
            }
        }
    };


private:

    // Private Data

        //- Maximum squared distance
        scalar maxDistance2_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const nearestFaceAMI&) = delete;

        autoPtr<mapDistribute> calcFaceMap
        (
            const List<nearestAndDist>& localInfo,
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch
        ) const;

        autoPtr<mapDistribute> calcDistributed
        (
            const primitivePatch& src,
            const primitivePatch& tgt,
            labelListList& srcToTgtAddr,
            scalarListList& srcToTgtWght
        ) const;


public:

    //- Runtime type information
    TypeName("nearestFaceAMI");


    // Constructors

        //- Construct from dictionary
        nearestFaceAMI
        (
            const dictionary& dict,
            const bool reverseTarget = false
        );

        //- Construct from components
        nearestFaceAMI
        (
            const bool requireMatch = true,
            const bool reverseTarget = false,
            const scalar lowWeightCorrection = -1
        );

        //- Construct as copy
        nearestFaceAMI(const nearestFaceAMI& ami);

        //- Construct and return a clone
        virtual autoPtr<AMIInterpolation> clone() const
        {
            return autoPtr<AMIInterpolation>(new nearestFaceAMI(*this));
        }


    //- Destructor
    virtual ~nearestFaceAMI() = default;


    // Member Functions

        //- Update addressing and weights
        virtual bool calculate
        (
            const primitivePatch& srcPatch,
            const primitivePatch& tgtPatch,
            const autoPtr<searchableSurface>& surfPtr = nullptr
        );

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

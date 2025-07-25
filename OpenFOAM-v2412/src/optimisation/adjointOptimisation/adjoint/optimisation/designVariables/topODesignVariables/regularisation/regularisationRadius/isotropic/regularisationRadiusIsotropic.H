/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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
    Foam::isotropic

Description
    An isotropic regularisationRadius (same in all spatial directions)


SourceFiles
    isotropic.C

\*---------------------------------------------------------------------------*/

#ifndef isotropic_H
#define isotropic_H

#include "regularisationRadius.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class isotropic Declaration
\*---------------------------------------------------------------------------*/

class isotropic
:
    public regularisationRadius
{

private:

    // Private Member Functions

        //- No copy construct
        isotropic(const isotropic&) = delete;

        //- No copy assignment
        void operator=(const isotropic&) = delete;


protected:

    // Protected data

        //- Multiplier of the wall thickness, used to obtain the radius for
        //- the second regularisation in bi-fluid topO
        scalar wallThicknessMult_;

        //- Smoothing radius of the first regulatisation
        dimensionedScalar radius_;


    // Protected Member Functions

        //- Compute smoothing radius, if not directly given
        scalar computeRadius(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("isotropic");


    // Constructors

        //- Construct from components
        isotropic
        (
            const fvMesh& mesh,
            const dictionary& dict,
            bool adjustWallThickness
        );


    //- Destructor
    virtual ~isotropic() = default;


    // Member Functions

        //- Add a Laplacian term with an isotropic diffusivity
        virtual void addRegularisationTerm
        (
            fvScalarMatrix& matrix,
            bool isTopoField
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

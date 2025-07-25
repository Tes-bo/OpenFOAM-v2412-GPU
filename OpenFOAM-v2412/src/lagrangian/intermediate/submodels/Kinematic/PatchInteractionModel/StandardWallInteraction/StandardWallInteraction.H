/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::StandardWallInteraction

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Wall interaction model.

    Three choices:
      - rebound - optionally specify elasticity and restitution coefficients
      - stick   - particles assigned zero velocity
      - escape  - remove particle from the domain

    Example usage:
    \verbatim
    StandardWallInteractionCoeffs
    {
        type        rebound; // stick, escape
        e           1;       // optional - elasticity coeff
        mu          0;       // optional - (1 - restitution coeff)
    }
    \endverbatim

Note
  - \c mu=1 results in perfectly inelastic collisions.

\*---------------------------------------------------------------------------*/

#ifndef StandardWallInteraction_H
#define StandardWallInteraction_H

#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class StandardWallInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class StandardWallInteraction
:
    public PatchInteractionModel<CloudType>
{
protected:

    // Protected data

        // Reference to mesh
        const fvMesh& mesh_;

        //- Interaction type
        typename PatchInteractionModel<CloudType>::interactionType
            interactionType_;

        //- Elasticity coefficient
        scalar e_;

        //- The unity minus the restitution coefficient
        scalar mu_;


        // Bookkeeping for particle fates

            //- Number of parcels escaped
            List<List<label>> nEscape_;

            //- Mass of parcels escaped
            List<List<scalar>> massEscape_;

            //- Number of parcels stuck to patches
            List<List<label>> nStick_;

            //- Mass of parcels stuck to patches
            List<List<scalar>> massStick_;


        //- Flag to output escaped/mass particles sorted by injectorID
        bool outputByInjectorId_;

        //- InjectorId to index map, when outputting escaped/stick/...
        //  particles sorted by injectorID
        Map<label> injIdToIndex_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("standardWallInteraction");


    // Constructors

        //- Construct from dictionary
        StandardWallInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        StandardWallInteraction(const StandardWallInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new StandardWallInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~StandardWallInteraction() = default;


    // Member Functions

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );


        //- Write patch interaction info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "StandardWallInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::RecycleInteraction

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Patch interaction model to perform
    parcel transfer and recycle from one patch to another.

    Parcels that hit a given 'outflow' patch
    are recycled to a given 'inflow' patch, with optional:
    - recycle fraction [0-1]

Usage
    Minimal example by using
    \c constant/reactingCloud1Properties.subModels.multiInteractionCoeffs
    \verbatim
    multiInteractionCoeffs
    {
        ...

        model1
        {
            // Mandatory entries (unmodifiable)
            patchInteractionModel    recycleInteraction;

            recycleInteractionCoeffs
            {
                recyclePatches
                (
                    (<outletPatch1> <inletPatch1>)
                    (<outletPatch2> <inletPatch2>)
                    ...
                );

                recycleFraction      0.8; // [0-1]

                // Optional entries (unmodifiable)
                outputByInjectorId   false;
            }
        }

        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                      | Type | Reqd | Dflt
      patchInteractionModel | Type name: recycleInteraction  <!--
                                                  --> | word | yes  | -
      recyclePatches | Names of outlet-inlet patch pairs     <!--
                                                -->   | (word word) | yes | -
      recycleFraction | Fraction of parcels recycled from <!--
                                       --> outlet to inlet | scalar | yes | -
      outputByInjectorId | Flag to output escaped/mass    <!--
                     --> particles sorted by injectorID | bool | no | false
    \endtable

SourceFiles
    RecycleInteraction.C

\*---------------------------------------------------------------------------*/

#ifndef RecycleInteraction_H
#define RecycleInteraction_H

#include "patchInjectionBase.H"
#include "IDLList.H"
#include "PtrList.H"
#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class RecycleInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class RecycleInteraction
:
    public PatchInteractionModel<CloudType>
{
protected:

    typedef typename CloudType::parcelType parcelType;

    // Protected Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Outlet-inlet patch pair to apply parcel recycling
        List<Pair<word>> recyclePatches_;

        //- Patch IDs of recyclePatches
        List<Pair<label>> recyclePatchesIds_;

        //- Parcel IDs of recycled parcels
        List<IDLList<parcelType>> recycledParcels_;


        // Bookkeeping for particle fates

            //- Number of parcels removed
            List<List<label>> nRemoved_;

            //- Mass of parcels removed
            List<List<scalar>> massRemoved_;

            //- Number of parcels injected
            List<List<label>> nInjected_;

            //- Mass of parcels injected
            List<List<scalar>> massInjected_;


        //- Injector ID to local index map
        Map<label> injIdToIndex_;

        //- Injection patch pointer
        PtrList<patchInjectionBase> injectionPatchPtr_;

        //- Parcel fraction to be recycled from outlet to inlet
        const scalar recycleFraction_;

        //- Flag to output escaped/mass particles sorted by injectorID
        bool outputByInjectorId_;


    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("recycleInteraction");


    // Constructors

        //- Construct from dictionary
        RecycleInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        RecycleInteraction(const RecycleInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new RecycleInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~RecycleInteraction() = default;


    // Member Functions

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );

        //- Post-evolve hook
        virtual void postEvolve();

        //- Write patch interaction info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RecycleInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

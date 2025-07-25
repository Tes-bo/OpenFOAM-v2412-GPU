/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd
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
    Foam::DTRMParticle

Description
    Discrete Transfer Radiation Model (DTRM) particle

SourceFiles
    DTRMParticle.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_DTRMParticle_H
#define Foam_DTRMParticle_H

#include "particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCell.H"
#include "volFieldsFwd.H"
#include "reflectionModel.H"
#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class DTRMParticle;

Ostream& operator<<(Ostream&, const DTRMParticle&);

using namespace Foam::radiation;

/*---------------------------------------------------------------------------*\
                        Class DTRMParticle Declaration
\*---------------------------------------------------------------------------*/

class DTRMParticle
:
    public particle
{
    // Private data

        //- Initial position
        point p0_;

        //- Target position
        point p1_;

        //- Initial radiation intensity [W/m2]
        scalar I0_;

        //- Radiation intensity [W/m2]
        scalar I_;

        //- Area of radiation
        scalar dA_;

        //- Trasnmissive index
        label transmissiveId_;


public:

    friend class Cloud<DTRMParticle>;

    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCell<scalar>& aInterp_;
            const interpolationCell<scalar>& eInterp_;
            const interpolationCell<scalar>& EInterp_;
            const interpolationCell<scalar>& TInterp_;
            const interpolationCellPoint<vector>& nHatInterp_;

            //- Reflected cells
            const labelField& relfectedCells_;

            //- Ptr to reflectiom model
            UPtrList<reflectionModel> reflection_;

            //- Heat source term
            volScalarField& Q_;


    public:

        // Constructors

            inline trackingData
            (
                Cloud<DTRMParticle>& spc,
                const interpolationCell<scalar>& aInterp,
                const interpolationCell<scalar>& eInterp,
                const interpolationCell<scalar>& EInterp,
                const interpolationCell<scalar>& TInterp,
                const interpolationCellPoint<vector>& nHatInterp,
                const labelField&,
                const UPtrList<reflectionModel>&,
                volScalarField& Q
            );


        // Member Functions

            inline const interpolationCell<scalar>& aInterp() const;
            inline const interpolationCell<scalar>& eInterp() const;
            inline const interpolationCell<scalar>& EInterp() const;
            inline const interpolationCell<scalar>& TInterp() const;
            inline const interpolationCellPoint<vector>& nHatInterp() const;
            inline const labelField& relfectedCells() const;
            inline const UPtrList<reflectionModel>& reflection() const;

            inline scalar& Q(label celli);
    };


    // Static Data Members

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


        //- String representation of properties
        AddToPropertyList
        (
            particle,
            " p0"
          + " p1"
          + " I0"
          + " I"
          + " dA"
          + " transmissiveId";
        );


    // Constructors

        //- Construct from components, with searching for tetFace and
        //  tetPt unless disabled by doCellFacePt = false.
        DTRMParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& targetPosition,
            const scalar I,
            const label cellI,
            const scalar dA,
            const label transmissiveId
        );

        //- Construct from components
        DTRMParticle
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti,
            const vector& position,
            const vector& targetPosition,
            const scalar I,
            const scalar dA,
            const label transmissiveId
        );

        //- Construct from Istream
        DTRMParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as copy
        DTRMParticle(const DTRMParticle& p);

        //- Return a clone
        virtual autoPtr<particle> clone() const
        {
            return particle::Clone(*this);
        }

        //- Factory class to read-construct particles (for parallel transfer)
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<DTRMParticle> operator()(Istream& is) const
            {
                return autoPtr<DTRMParticle>
                (
                    new DTRMParticle(mesh_, is, true)
                );
            }
        };


    // Access

        //- Return const access to the initial position
        const point& p0() const noexcept { return p0_; }

        //- Return const access to the target position
        const point& p1() const noexcept { return p1_; }

        //- Return const access to the initial intensity
        scalar I0() const noexcept { return I0_; }

        //- Return const access to the current intensity
        scalar I() const noexcept { return I_; }

        //- Return const access dA
        scalar dA() const noexcept { return dA_; }


    // Edit

        //- Return access to the target position
        point& p1() noexcept { return p1_; }

        //- Return access to the initial intensity
        scalar& I0() noexcept { return I0_; }

        //- Return access to the current intensity
        scalar& I() noexcept { return I_; }

        //- Return access to dA
        scalar& dA() noexcept { return dA_; }


   // Tracking

        //- Move
        bool move(Cloud<DTRMParticle>& , trackingData&, const scalar);


    // Member Operators

        //- Overridable function to handle the particle hitting a processorPatch
        void hitProcessorPatch
        (
            Cloud<DTRMParticle>&,
            trackingData& td
        );

        //- Overridable function to handle the particle hitting a wallPatch
        void hitWallPatch
        (
            Cloud<DTRMParticle>&,
            trackingData& td
        );

        bool hitPatch
        (
            Cloud<DTRMParticle>&,
            trackingData& td
        );


        // I-O

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly = false
            ) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const DTRMParticle& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DTRMParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

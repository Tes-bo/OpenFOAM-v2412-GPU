/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Norbert Weber, HZDR
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
    Foam::weightedFlux

Description
    Weighted flux interpolation scheme class.

    This scheme is used to compute fluxes with variable diffusivity or
    conductivity, as e.g.
    - a thermal flux: lambda*grad(T)
    - a mass flux: D*grad(u)
    - an electric current: -sigma*grad(potential)

    When using the Gauss theorem to compute a gradient, cell centred values
    need to be interpolated to the faces. Using this scheme, temperature (T)
    is weighted by thermal conductivity when being interpolated. Similarly,
    velocity is weighted by diffusivity (D) and the electric potential by
    the electric conductivity (sigma). Lambda, D or sigma are read from the
    object registry - the names need to be specified in fvSchemes as e.g.

    \verbatim
        gradSchemes
        {
            grad(T)             Gauss weightedFlux lambda;
            grad(u)             Gauss weightedFlux D;
            grad(potential)     Gauss weightedFlux sigma;
        }
    \endverbatim

    For more details, see equation 16 and 17 in
    \verbatim
        Weber, N., Beckstein, P., Galindo, V., Starace, M. & Weier, T. (2018).
        Electro-vortex flow simulation using coupled meshes.
        Computers and Fluids 168, 101-109.
        doi:10.1016/j.compfluid.2018.03.047
        https://arxiv.org/pdf/1707.06546.pdf
    \endverbatim

Note
    For support, contact Norbert.Weber@hzdr.de

SourceFiles
    weightedFlux.C

\*---------------------------------------------------------------------------*/

#ifndef weightedFlux_H
#define weightedFlux_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class weightedFlux Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class weightedFlux
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Const reference to step-wise pre-gradient factor field
        const volScalarField& sigma_;

    // Demand-driven data

        //- Face to owner cell distance
        mutable surfaceScalarField* oDelta_;

        //- Face to neighbour cell distance
        mutable surfaceScalarField* nDelta_;


    // Private Member Functions

        //- Compute face-owner and face-neighbour distance
        void makeDeltas() const;

        //- No copy assignment
        void operator=(const weightedFlux&) = delete;


protected:

    // Protected Member Functions

        // Storage management

            //- Clear all fields
            void clearOut();


public:

    //- Runtime type information
    TypeName("weightedFlux");


    // Constructors

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        weightedFlux
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            sigma_(this->mesh().objectRegistry::template
                lookupObject<volScalarField>(word(is))),
            oDelta_(nullptr),
            nDelta_(nullptr)
        {}

        //- Construct from faceFlux and Istream
        weightedFlux
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            sigma_(this->mesh().objectRegistry::template
                lookupObject<volScalarField>(word(is))),
            oDelta_(nullptr),
            nDelta_(nullptr)
        {}


    //- Destructor
    virtual ~weightedFlux();


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().surfaceInterpolation::weights();
        }

        //- Return the distance between face and owner cell
        const surfaceScalarField& oDelta() const
        {
            if (!oDelta_)
            {
                makeDeltas();
            }

            return *oDelta_;
        }

        //- Return the distance between face and neighbour cell
        const surfaceScalarField& nDelta() const
        {
            if (!nDelta_)
            {
                makeDeltas();
            }

            return *nDelta_;
        }

        //- Interpolate the cell values to faces
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::sensitivitySurface

Description
    Calculation of adjoint based sensitivities at wall faces

    Reference:
    \verbatim
        The computation of the sensitivity derivatives  follows the (E-)SI
        formulation of

            Kavvadias, I. S., Papoutsis-Kiachagias, E. M.,
            Giannakoglou, K. C. (2015).
            On the proper treatment of grid sensitivities in continuous adjoint
            methods for shape optimization.
            Journal of computational physics, 301, 1-18.
            https://doi.org/10.1016/j.jcp.2015.08.012

        whereas their smoothing based on the computation of the 'Sobolev
        gradient' is derived from

            Vassberg J. C., Jameson A. (2006).
            Aerodynamic Shape Optimization Part I: Theoretical Background.
            VKI Lecture Series, Introduction to Optimization and
            Multidisciplinary Design, Brussels, Belgium, 8 March, 2006.
    \endverbatim

SourceFiles
    sensitivitySurface.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivitySurfaceIncompressible_H
#define sensitivitySurfaceIncompressible_H

#include "sensitivitySurfacePoints.H"
#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class sensitivitySurface Declaration
\*---------------------------------------------------------------------------*/

class sensitivitySurface
:
    public sensitivitySurfacePoints
{
protected:

    // Protected data

        //- Smooth sensitivity derivatives based on a surface Laplace solver
        bool smoothSensitivities_;

        //- Return the complete vector of sensitivities
        bool returnVectorField_;

        //bool finalResultIncludesArea_;


    // Protected Member Functions

        //- Compute the size of the return field
        label computeFaceDerivativesSize(const bool computeVectorFieldSize);

        //- Smooth sensitivity derivatives based on the computation of the
        //- 'Sobolev gradient'
        void smoothSensitivities();

        //- Compute the physical smoothing radius based on the average boundary
        //- face 'length'
        scalar computeRadius(const faMesh& aMesh);


private:

    // Private Member Functions

        //- No copy construct
        sensitivitySurface(const sensitivitySurface&) = delete;

        //- No copy assignment
        void operator=(const sensitivitySurface&) = delete;


public:

    //- Runtime type information
    TypeName("surface");


    // Constructors

        //- Construct from components
        sensitivitySurface
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~sensitivitySurface() = default;


    // Member Functions

        //- Read controls and update solver pointers if necessary
        void read();

        //- Assemble sensitivities
        virtual void assembleSensitivities
        (
            autoPtr<designVariables>& designVars
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

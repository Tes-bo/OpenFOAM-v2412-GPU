/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::conjugateGradient

Description
    The Conjugate Gradient formula.

    Beta computed from Fletcher-Reeves or Polak-Ribiere approaches.
    Polak-Ribiere seems to misbehave.
    Not tested with line search.

SourceFiles
    conjugateGradient.C

\*---------------------------------------------------------------------------*/

#ifndef conjugateGradient_H
#define conjugateGradient_H

#include "updateMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class conjugateGradient Declaration
\*---------------------------------------------------------------------------*/

class conjugateGradient
:
    public updateMethod
{
protected:

    // Protected data

        scalarField dxOld_;
        scalarField sOld_;
        word betaType_;


private:

    // Private Member Functions

        //- No copy construct
        conjugateGradient(const conjugateGradient&) = delete;

        //- No copy assignment
        void operator=(const conjugateGradient&) = delete;


public:

    //- Runtime type information
    TypeName("conjugateGradient");


    // Constructors

        //- Construct from components
        conjugateGradient
        (
            const fvMesh& mesh,
            const dictionary& dict,
            autoPtr<designVariables>& designVars,
            const label nConstraints,
            const word& type
        );


    //- Destructor
    virtual ~conjugateGradient() = default;


    // Member Functions

        //- Compute design variables correction
        void computeCorrection();

        //- Update old correction. For use when eta has been changed externally
        virtual void updateOldCorrection(const scalarField& oldCorrection);

        //- Write useful quantities to files
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

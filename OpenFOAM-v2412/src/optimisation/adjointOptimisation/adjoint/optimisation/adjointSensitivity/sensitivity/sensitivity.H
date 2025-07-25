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
    Foam::sensitivity

Description
    Abstract base class for adjoint sensitivities

SourceFiles
    sensitivity.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivity_H
#define sensitivity_H

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "volPointInterpolation.H"

#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointPatchFieldsFwd.H"
#include "fixedValuePointPatchField.H"
#include "boundaryFieldsFwd.H"
#include "createZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class designVariables;

/*---------------------------------------------------------------------------*\
                         Class sensitivity Declaration
\*---------------------------------------------------------------------------*/

class sensitivity
{
protected:

    // Protected data

        const fvMesh& mesh_;
        dictionary dict_;
        bool writeFieldSens_;

        // Field sensitivities. Topology optimisation
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<volScalarField> fieldSensPtr_;


private:

    // Private Member Functions

        //- No copy construct
        sensitivity(const sensitivity&) = delete;

        //- No copy assignment
        void operator=(const sensitivity&) = delete;


public:

    //- Runtime type information
    TypeName("sensitivity");

    // Constructors

        //- Construct from components
        sensitivity
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~sensitivity() = default;


    // Member Functions

        //- Return reference to mesh
        inline const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the construction dictionary
        inline const dictionary& dict() const
        {
            return
                dict_.optionalSubDict(mesh_.name()).
                    optionalSubDict("sensitivities");
        }

        //- Read dictionary if changed
        virtual bool readDict(const dictionary& dict);

        //- Calculates and returns sensitivity field
        virtual const scalarField& calculateSensitivities
        (
            autoPtr<designVariables>& designVars
        ) = 0;

        //- Get the fieldSensPtr
        inline const autoPtr<volScalarField>& fieldSensPtr() const
        {
            return fieldSensPtr_;
        }

        //- Write sensitivity fields.
        //  If valid, copies boundaryFields to volFields and writes them.
        //  Virtual to be reimplemented by control points-based methods
        //  (Bezier, RBF) which do not need to write fields
        virtual void write(const word& baseName = word::null);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "sensitivityTemplates.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

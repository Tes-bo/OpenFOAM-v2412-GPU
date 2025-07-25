/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::functionObjects::energySpectrum

Group
    grpFieldFunctionObjects

Description
    Calculates the energy spectrum for a structured IJK mesh

Usage
    Example of function object specification:
    \verbatim
    energySpectrum1
    {
        type        energySpectrum;
        libs        (fieldFunctionObjects);
        ...
        log         yes;
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description               | Required    | Default value
        type         | type name: energySpectrum | yes         |
        log          | write info to standard output | no      | yes
    \endtable

    Output data is written to the file \<timeDir\>/energySpectrum.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    energySpectrum.C
    energySpectrumTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_energySpectrum_H
#define functionObjects_energySpectrum_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "Vector.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class energySpectrum Declaration
\*---------------------------------------------------------------------------*/

class energySpectrum
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- I-J-K mesh addressing
        labelList cellAddr_;

        //- Name of velocity field, default = U
        word UName_;

        //- Number of cells in I-J-K directions
        Vector<int> N_;

        //- Reference point
        vector c0_;

        //- Cell length scale
        vector deltaC_;

        //- Wave number
        scalar kappaNorm_;


    // Protected Member Functions


        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Calculate and write the spectrum
        void calcAndWriteSpectrum
        (
            const vectorField& U,
            const vectorField& C,
            const vector& c0,
            const vector& deltaC,
            const Vector<int>& N,
            const scalar kappaNorm
        );

        //- No copy construct
        energySpectrum(const energySpectrum&) = delete;

        //- No copy assignment
        void operator=(const energySpectrum&) = delete;


public:

    //- Runtime type information
    TypeName("energySpectrum");


    // Constructors

        //- Construct from Time and dictionary
        energySpectrum
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~energySpectrum() = default;


    // Member Functions

        //- Read the field min/max data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the energySpectrum
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::windowModel

Description
    Base class for windowing models

SourceFiles
    windowModel.C

\*---------------------------------------------------------------------------*/

#ifndef windowModel_H
#define windowModel_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class windowModel Declaration
\*---------------------------------------------------------------------------*/

class windowModel
:
    public scalarField
{

protected:

    // Protected Data

        //- Overlap percent
        scalar overlapPercent_;

        //- Number of overlap samples per window
        label nOverlapSamples_;

        //- Number of windows
        label nWindow_;


public:

    //- Runtime type information
    TypeName("windowModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        windowModel,
        dictionary,
        (
            const dictionary& dict,
            const label nSamples
        ),
        (dict, nSamples)
    );


    //- Construct from dictionary
    windowModel(const dictionary& dict, const label nSamples);


    // Selectors

        //- Return a reference to the selected window model
        static autoPtr<windowModel> New
        (
            const dictionary& dict,
            const label nSamples
        );


    //- Destructor
    virtual ~windowModel() = default;


    // Public Member Functions

        //- Return the number of samples in the window
        label nSamples() const;

        //- Return the overlap percent
        scalar overlapPercent() const;

        //- Return number of overlap samples per window
        label nOverlapSamples() const;

        //- Return the number of windows
        label nWindow() const;

        //- Return the total number of windows for a given number of samples
        label nWindowsTotal(label nSamplesTotal) const;

        //- Validate that the window is applicable to the data set size, and
        //  return the number of required data points
        label validate(label n);

        //- Return the windowed data
        template<class Type>
        tmp<Field<Type>> apply
        (
            const Field<Type>& fld,
            const label windowI
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "windowModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

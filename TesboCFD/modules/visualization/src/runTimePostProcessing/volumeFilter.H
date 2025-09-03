/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::runTimePostPro::volumeFilter

Description
    Visualisation of OpenFOAM volume fields as surface data using
    a VTK filter cascade.

Note
    Since this filter includes an OpenFOAM/VTK adaptor level,
    it is ill-suited to mismatches in data parallelization.
    If OpenFOAM is running in parallel but VTK is not, it would be rather
    expensive to collect all the data on the master node for this filter.
    That approach is acceptable for smaller amounts of data, but do not
    allow for volume meshes.

SourceFiles
    volumeFilter.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_volumeFilter_H
#define functionObjects_runTimePostPro_volumeFilter_H

#include "runTimePostProcessing.H"
#include "surface.H"
#include "foamVtkVtuAdaptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkDataSet;
class vtkMultiPieceDataSet;
template<class T> class vtkSmartPointer;

namespace Foam
{
namespace vtk
{
class vtuAdaptor;
}

namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                        Class volumeFilter Declaration
\*---------------------------------------------------------------------------*/

class volumeFilter
:
    public surface
{
protected:

    // Protected Member Functions

        //- Return a vtu mesh with addressing information stored in adaptor
        vtkSmartPointer<vtkMultiPieceDataSet> mesh
        (
            Foam::vtk::vtuAdaptor& adaptor
        ) const;


        bool addDimField
        (
            vtkDataSet* piece,
            const vtk::vtuAdaptor& adaptor,
            const regIOobject* ioptr,
            const word& fieldName
        ) const;

        int addDimField
        (
            vtkMultiPieceDataSet* multiPiece,
            const vtk::vtuAdaptor& adaptor,
            const regIOobject* ioptr,
            const word& fieldName
        ) const;

        template<class Type>
        bool addDimField
        (
            vtkDataSet* piece,
            const vtk::vtuAdaptor& adaptor,
            const regIOobject* ioptr,
            const word& fieldName
        ) const;

        template<class Type>
        int addDimField
        (
            vtkMultiPieceDataSet* multiPiece,
            const vtk::vtuAdaptor& adaptor,
            const regIOobject* ioptr,
            const word& fieldName
        ) const;


        //- No copy construct
        volumeFilter(const volumeFilter&) = delete;

        //- No copy assignment
        void operator=(const volumeFilter&) = delete;


public:

    // Constructors

        //- Construct from dictionary
        volumeFilter
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~volumeFilter() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volumeFilterTemplates.C"
#endif

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    temporalInterpolate

Group
    grpPostProcessingUtilities

Description
    Interpolate fields between time-steps e.g. for animation.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"

#include "fvMesh.H"
#include "Time.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "ReadFields.H"
#include "interpolationWeights.H"
#include "uniformInterpolate.H"

using namespace Foam;

class fieldInterpolator
{
    Time& runTime_;
    const fvMesh& mesh_;
    const IOobjectList& objects_;
    const wordRes& selectedFields_;
    instant ti_;
    instant ti1_;
    const interpolationWeights& interpolator_;
    const wordList& timeNames_;
    int divisions_;

public:

    fieldInterpolator
    (
        Time& runTime,
        const fvMesh& mesh,
        const IOobjectList& objects,
        const wordRes& selectedFields,
        const instant& ti,
        const instant& ti1,
        const interpolationWeights& interpolator,
        const wordList& timeNames,
        int divisions
    )
    :
        runTime_(runTime),
        mesh_(mesh),
        objects_(objects),
        selectedFields_(selectedFields),
        ti_(ti),
        ti1_(ti1),
        interpolator_(interpolator),
        timeNames_(timeNames),
        divisions_(divisions)
    {}

    template<class GeoFieldType>
    void interpolate();
};


template<class GeoFieldType>
void fieldInterpolator::interpolate()
{
    label nFields = 0;
    for
    (
        const IOobject& io :
        (
            selectedFields_.empty()
          ? objects_.csorted<GeoFieldType>()
          : objects_.csorted<GeoFieldType>(selectedFields_)
        )
    )
    {
        const word& fieldName = io.name();

        if (!nFields)
        {
            Info<< "    " << GeoFieldType::typeName << 's';
        }

        Info<< ' ' << fieldName << '(';

        const scalar deltaT = (ti1_.value() - ti_.value())/(divisions_ + 1);

        for (int j=0; j<divisions_; j++)
        {
            instant timej = instant(ti_.value() + (j + 1)*deltaT);

            runTime_.setTime(instant(timej.name()), 0);

            Info<< timej.name();

            if (j < divisions_-1)
            {
                Info<< " ";
            }

            // Calculate times to read and weights
            labelList indices;
            scalarField weights;
            interpolator_.valueWeights
            (
                runTime_.value(),
                indices,
                weights
            );

            const wordList selectedTimeNames
            (
                UIndirectList<word>(timeNames_, indices)()
            );

            //Info<< "For time " << runTime_.value()
            //    << " need times " << selectedTimeNames
            //    << " need weights " << weights << endl;


            // Read on the objectRegistry all the required fields
            ReadFields<GeoFieldType>
            (
                fieldName,
                mesh_,
                selectedTimeNames
            );

            GeoFieldType fieldj
            (
                uniformInterpolate<GeoFieldType>
                (
                    IOobject
                    (
                        fieldName,
                        runTime_.timeName(),
                        io.db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        IOobject::NO_REGISTER
                    ),
                    fieldName,
                    selectedTimeNames,
                    weights
                )
            );

            fieldj.write();
        }

        Info<< ')';
        ++nFields;
    }

    if (nFields) Info<< endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Interpolate fields between time-steps. Eg, for animation."
    );

    timeSelector::addOptions();
    #include "addRegionOption.H"
    argList::addOption
    (
        "fields",
        "wordRes",
        "The fields (or field) to be interpolated."
        " Eg, '(U T p \"Y.*\")' or a single field 'U'"
    );
    argList::addOption
    (
        "divisions",
        "integer",
        "Specify number of temporal sub-divisions to create (default = 1)."
    );
    argList::addOption
    (
        "interpolationType",
        "word",
        "The type of interpolation (linear or spline)"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"

    // Non-mandatory
    const wordRes selectedFields(args.getList<wordRe>("fields", false));

    if (selectedFields.empty())
    {
        Info<< "Interpolating all fields" << nl << endl;
    }
    else
    {
        Info<< "Interpolating fields " << flatOutput(selectedFields)
            << nl << endl;
    }


    const int divisions = args.getOrDefault<int>("divisions", 1);
    Info<< "Using " << divisions << " per time interval" << nl << endl;


    const word interpolationType =
        args.getOrDefault<word>("interpolationType", "linear");

    Info<< "Using interpolation " << interpolationType << nl << endl;


    instantList timeDirs = timeSelector::select0(runTime, args);

    scalarField timeVals(timeDirs.size());
    wordList timeNames(timeDirs.size());
    forAll(timeDirs, i)
    {
        timeVals[i] = timeDirs[i].value();
        timeNames[i] = timeDirs[i].name();
    }
    autoPtr<interpolationWeights> interpolatorPtr
    (
        interpolationWeights::New
        (
            interpolationType,
            timeVals
        )
    );


    #include "createNamedMesh.H"

    Info<< "Interpolating fields for times:" << endl;

    for (label timei = 0; timei < timeDirs.size() - 1; timei++)
    {
        runTime.setTime(timeDirs[timei], timei);

        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());

        fieldInterpolator interpolator
        (
            runTime,
            mesh,
            objects,
            selectedFields,
            timeDirs[timei],
            timeDirs[timei+1],
            interpolatorPtr(),
            timeNames,
            divisions
        );

        // Interpolate vol fields
        interpolator.interpolate<volScalarField>();
        interpolator.interpolate<volVectorField>();
        interpolator.interpolate<volSphericalTensorField>();
        interpolator.interpolate<volSymmTensorField>();
        interpolator.interpolate<volTensorField>();

        // Interpolate surface fields
        interpolator.interpolate<surfaceScalarField>();
        interpolator.interpolate<surfaceVectorField>();
        interpolator.interpolate<surfaceSphericalTensorField>();
        interpolator.interpolate<surfaceSymmTensorField>();
        interpolator.interpolate<surfaceTensorField>();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

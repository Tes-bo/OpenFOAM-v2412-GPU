/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

Info << nl << "Creating autocorrelation functions." << endl;

IOdictionary mdTransportProperitesDict
(
    IOobject
    (
        "mdTransportProperitesDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    )
);

const dictionary& autocorrelationFunctionDict
(
    mdTransportProperitesDict.subDict("autocorrelationFunctions")
);

//- Velocity autocorrelation function

Info << tab << "velocty" << endl;

const dictionary& velocityACFDict
(
    autocorrelationFunctionDict.subDict("velocity")
);

correlationFunction<vector> vacf
(
    mesh,
    velocityACFDict,
    molecules.size()
);

bool writeVacf(velocityACFDict.get<bool>("writeFile"));

//- Pressure autocorrelation function

Info << tab << "pressure" << endl;

const dictionary& pressureACFDict
(
    autocorrelationFunctionDict.subDict("pressure")
);

correlationFunction<vector> pacf
(
    mesh,
    pressureACFDict,
    1
);

bool writePacf(pressureACFDict.get<bool>("writeFile"));

//- Heat flux autocorrelation function

Info << tab << "heat flux" << endl;

const dictionary& heatFluxACFDict
(
    autocorrelationFunctionDict.subDict("heatFlux")
);

correlationFunction<vector> hfacf
(
    mesh,
    heatFluxACFDict,
    1
);

bool writeHFacf(heatFluxACFDict.get<bool>("writeFile"));

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "layerParameters.H"
#include "polyBoundaryMesh.H"
#include "unitConversion.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "medialAxisMeshMover.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::layerParameters::thicknessModelType
>
Foam::layerParameters::thicknessModelTypeNames_
({
    { thicknessModelType::FIRST_AND_TOTAL, "firstAndOverall" },
    { thicknessModelType::FIRST_AND_EXPANSION, "firstAndExpansion" },
    { thicknessModelType::FINAL_AND_TOTAL, "finalAndOverall" },
    { thicknessModelType::FINAL_AND_EXPANSION, "finalAndExpansion" },
    { thicknessModelType::TOTAL_AND_EXPANSION, "overallAndExpansion" },
    { thicknessModelType::FIRST_AND_RELATIVE_FINAL, "firstAndRelativeFinal" },
});

const Foam::scalar Foam::layerParameters::defaultConcaveAngle = 90;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::layerParameters::layerExpansionRatio
(
    const label n,
    const scalar totalOverFirst
)
{
    if (n <= 1)
    {
        return 1.0;
    }

    const label maxIters = 20;
    const scalar tol = 1e-8;

    if (mag(n-totalOverFirst) < tol)
    {
        return 1.0;
    }

    // Calculate the bounds of the solution
    scalar minR;
    scalar maxR;

    if (totalOverFirst < n)
    {
        minR = 0.0;
        maxR = pow(totalOverFirst/n, scalar(1)/(n-1));
    }
    else
    {
        minR = pow(totalOverFirst/n, scalar(1)/(n-1));
        maxR = totalOverFirst/(n - 1);
    }

    // Starting guess
    scalar r = 0.5*(minR + maxR);

    for (label i = 0; i < maxIters; ++i)
    {
        const scalar prevr = r;

        const scalar fx = pow(r, n) - totalOverFirst*r - (1 - totalOverFirst);
        const scalar dfx = n*pow(r, n - 1) - totalOverFirst;
        r -= fx/dfx;

        if (mag(r - prevr) < tol)
        {
            break;
        }
    }
    return r;
}


void Foam::layerParameters::readLayerParameters
(
    const bool verbose,
    const dictionary& dict,
    const thicknessModelType& spec,
    scalar& firstLayerThickness,
    scalar& finalLayerThickness,
    scalar& thickness,
    scalar& expansionRatio
)
{
    // Now we have determined the layer-specification read the actual fields
    switch (spec)
    {
        case FIRST_AND_TOTAL:
            if (verbose)
            {
                Info<< "Layer specification as" << nl
                    << "- first layer thickness ('firstLayerThickness')" << nl
                    << "- overall thickness ('thickness')" << endl;
            }
            firstLayerThickness = dict.get<scalar>("firstLayerThickness");
            thickness = dict.get<scalar>("thickness");
        break;

        case FIRST_AND_EXPANSION:
            if (verbose)
            {
                Info<< "Layer specification as" << nl
                    << "- first layer thickness ('firstLayerThickness')" << nl
                    << "- expansion ratio ('expansionRatio')" << endl;
            }
            firstLayerThickness = dict.get<scalar>("firstLayerThickness");
            expansionRatio = dict.get<scalar>("expansionRatio");
        break;

        case FINAL_AND_TOTAL:
            if (verbose)
            {
                Info<< "Layer specification as" << nl
                    << "- final layer thickness ('finalLayerThickness')" << nl
                    << "- overall thickness ('thickness')" << endl;
            }
            finalLayerThickness = dict.get<scalar>("finalLayerThickness");
            thickness = dict.get<scalar>("thickness");
        break;

        case FINAL_AND_EXPANSION:
            if (verbose)
            {
                Info<< "Layer specification as" << nl
                    << "- final layer thickness ('finalLayerThickness')" << nl
                    << "- expansion ratio ('expansionRatio')" << endl;
            }
            finalLayerThickness = dict.get<scalar>("finalLayerThickness");
            expansionRatio = dict.get<scalar>("expansionRatio");
        break;

        case TOTAL_AND_EXPANSION:
            if (verbose)
            {
                Info<< "Layer specification as" << nl
                    << "- overall thickness ('thickness')" << nl
                    << "- expansion ratio ('expansionRatio')" << endl;
            }
            thickness = dict.get<scalar>("thickness");
            expansionRatio = dict.get<scalar>("expansionRatio");
        break;

        case FIRST_AND_RELATIVE_FINAL:
            if (verbose)
            {
                Info<< "Layer specification as" << nl
                    << "- absolute first layer thickness"
                    << " ('firstLayerThickness')"
                    << nl
                    << "- and final layer thickness"
                    << " ('finalLayerThickness')" << nl
                    << endl;
            }
            firstLayerThickness = dict.get<scalar>("firstLayerThickness");
            finalLayerThickness = dict.get<scalar>("finalLayerThickness");
        break;

        default:
            FatalIOErrorInFunction(dict)
                << "problem." << exit(FatalIOError);
        break;
    }
}


void Foam::layerParameters::calculateLayerParameters
(
    const thicknessModelType& spec,
    const label nLayers,
    scalar& firstThickness,
    scalar& finalThickness,
    scalar& thickness,
    scalar& expansionRatio
)
{
    // Calculate the non-read parameters
    switch (spec)
    {
        case FIRST_AND_TOTAL:
            expansionRatio = layerExpansionRatio
            (
                spec,
                nLayers,
                firstThickness,
                VGREAT,
                thickness,              //totalThickness
                VGREAT                  //expansionRatio
            );
            finalThickness =
                thickness
               *finalLayerThicknessRatio(nLayers, expansionRatio);

        break;

        case FIRST_AND_EXPANSION:
            thickness = layerThickness
            (
                spec,
                nLayers,
                firstThickness,         //firstThickness
                VGREAT,                 //finalThickness
                VGREAT,                 //totalThickness
                expansionRatio          //expansionRatio
            );
            finalThickness =
                thickness
               *finalLayerThicknessRatio(nLayers, expansionRatio);

        break;

        case FINAL_AND_TOTAL:
            firstThickness = firstLayerThickness
            (
                spec,
                nLayers,
                VGREAT,                 //firstThickness
                VGREAT,                 //finalThickness
                thickness,              //totalThickness
                VGREAT                  //expansionRatio
            );
            expansionRatio = layerExpansionRatio
            (
                spec,
                nLayers,
                VGREAT,                 //firstThickness
                finalThickness,         //finalThickness
                thickness,              //totalThickness
                VGREAT                  //expansionRatio
            );

        break;

        case FINAL_AND_EXPANSION:
            firstThickness = firstLayerThickness
            (
                spec,
                nLayers,
                VGREAT,                 //firstThickness
                finalThickness,         //finalThickness
                VGREAT,                 //thickness
                expansionRatio          //expansionRatio
            );
            thickness = layerThickness
            (
                spec,
                nLayers,
                VGREAT,                 //firstThickness
                finalThickness,         //finalThickness
                VGREAT,                 //totalThickness
                expansionRatio          //expansionRatio
            );
        break;

        case TOTAL_AND_EXPANSION:
            firstThickness = firstLayerThickness
            (
                spec,
                nLayers,
                VGREAT,                 //firstThickness
                finalThickness,         //finalThickness
                VGREAT,                 //thickness
                expansionRatio          //expansionRatio
            );
            finalThickness =
                thickness
               *finalLayerThicknessRatio(nLayers, expansionRatio);

        break;

        case FIRST_AND_RELATIVE_FINAL:
            thickness = layerThickness
            (
                spec,
                nLayers,
                firstThickness,         //firstThickness
                finalThickness,         //finalThickness
                VGREAT,                 //totalThickness
                VGREAT                  //expansionRatio
            );
            expansionRatio = layerExpansionRatio
            (
                spec,
                nLayers,
                firstThickness,         //firstThickness
                finalThickness,         //finalThickness
                VGREAT,                 //totalThickness
                VGREAT                  //expansionRatio
            );

        break;

        default:
            FatalErrorInFunction << "Illegal thicknessModel " << spec
                << exit(FatalError);
        break;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh,
    const bool dryRun
)
:
    dict_(dict),
    dryRun_(dryRun),
    numLayers_(boundaryMesh.size(), -1),
    relativeSizes_
    (
        boundaryMesh.size(),
        meshRefinement::get<bool>(dict, "relativeSizes", dryRun)
    ),
    layerModels_(boundaryMesh.size(), FIRST_AND_TOTAL),
    firstLayerThickness_(boundaryMesh.size(), -123),
    finalLayerThickness_(boundaryMesh.size(), -123),
    thickness_(boundaryMesh.size(), -123),
    expansionRatio_(boundaryMesh.size(), -123),
    minThickness_
    (
        boundaryMesh.size(),
        meshRefinement::get<scalar>(dict, "minThickness", dryRun)
    ),
    featureAngle_(meshRefinement::get<scalar>(dict, "featureAngle", dryRun)),
    mergePatchFacesAngle_
    (
        dict.getOrDefault<scalar>
        (
            "mergePatchFacesAngle",
            featureAngle_
        )
    ),
    concaveAngle_
    (
        dict.getOrDefault<scalar>("concaveAngle", defaultConcaveAngle)
    ),
    nGrow_(meshRefinement::get<label>(dict, "nGrow", dryRun)),
    maxFaceThicknessRatio_
    (
        meshRefinement::get<scalar>(dict, "maxFaceThicknessRatio", dryRun)
    ),
    nBufferCellsNoExtrude_
    (
        meshRefinement::get<label>(dict, "nBufferCellsNoExtrude", dryRun)
    ),
    nLayerIter_(meshRefinement::get<label>(dict, "nLayerIter", dryRun)),
    nRelaxedIter_(labelMax),
    additionalReporting_(dict.getOrDefault("additionalReporting", false)),
    meshShrinker_
    (
        dict.getOrDefault
        (
            "meshShrinker",
            medialAxisMeshMover::typeName
        )
    ),
    nOuterIter_(dict.getOrDefault<scalar>("nOuterIter", 1))
{
    // Detect layer specification mode

    word spec;
    if (dict.readIfPresent("thicknessModel", spec))
    {
        layerModels_ = thicknessModelTypeNames_[spec];
    }
    else
    {
        // Count number of specifications
        label nSpec = 0;

        bool haveFirst = dict.found("firstLayerThickness");
        if (haveFirst)
        {
            nSpec++;
        }
        bool haveFinal = dict.found("finalLayerThickness");
        if (haveFinal)
        {
            nSpec++;
        }
        bool haveTotal = dict.found("thickness");
        if (haveTotal)
        {
            nSpec++;
        }
        bool haveExp = dict.found("expansionRatio");
        if (haveExp)
        {
            nSpec++;
        }

        if (nSpec == 2 && haveFirst && haveTotal)
        {
            layerModels_ = FIRST_AND_TOTAL;
            //Info<< "Layer thickness specified as first layer"
            //    << " and overall thickness." << endl;
        }
        else if (nSpec == 2 && haveFirst && haveExp)
        {
            layerModels_ = FIRST_AND_EXPANSION;
            //Info<< "Layer thickness specified as first layer"
            //    << " and expansion ratio." << endl;
        }
        else if (nSpec == 2 && haveFinal && haveTotal)
        {
            layerModels_ = FINAL_AND_TOTAL;
            //Info<< "Layer thickness specified as final layer"
            //    << " and overall thickness." << endl;
        }
        else if (nSpec == 2 && haveFinal && haveExp)
        {
            layerModels_ = FINAL_AND_EXPANSION;
            //Info<< "Layer thickness specified as final layer"
            //    << " and expansion ratio." << endl;
        }
        else if (nSpec == 2 && haveTotal && haveExp)
        {
            layerModels_ = TOTAL_AND_EXPANSION;
            //Info<< "Layer thickness specified as overall thickness"
            //    << " and expansion ratio." << endl;
        }
        else if (nSpec == 2 && haveFirst && haveFinal)
        {
            layerModels_ = FIRST_AND_RELATIVE_FINAL;
            //Info<< "Layer thickness specified as absolute first and"
            //    << " relative final layer ratio." << endl;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Over- or underspecified layer thickness."
                << " Please specify" << nl
                << "    first layer thickness ('firstLayerThickness')"
                << " and overall thickness ('thickness') or" << nl
                << "    first layer thickness ('firstLayerThickness')"
                << " and expansion ratio ('expansionRatio') or" << nl
                << "    final layer thickness ('finalLayerThickness')"
                << " and expansion ratio ('expansionRatio') or" << nl
                << "    final layer thickness ('finalLayerThickness')"
                << " and overall thickness ('thickness') or" << nl
                << "    overall thickness ('thickness')"
                << " and expansion ratio ('expansionRatio'"
                << exit(FatalIOError);
        }
    }


    // Now we have determined the layer-specification read the actual fields
    scalar firstThickness;
    scalar finalThickness;
    scalar thickness;
    scalar expansionRatio;
    readLayerParameters
    (
        true,               // verbose
        dict,
        layerModels_[0],    // spec
        firstThickness,
        finalThickness,
        thickness,
        expansionRatio
    );
    firstLayerThickness_ = firstThickness;
    finalLayerThickness_ = finalThickness;
    thickness_ = thickness;
    expansionRatio_ = expansionRatio;

    dict.readIfPresent("nRelaxedIter", nRelaxedIter_);

    if (nLayerIter_ < 0 || nRelaxedIter_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Layer iterations should be >= 0." << nl
            << "nLayerIter:" << nLayerIter_
            << " nRelaxedIter:" << nRelaxedIter_
            << exit(FatalIOError);
    }


    const dictionary& layersDict = meshRefinement::subDict
    (
        dict,
        "layers",
        dryRun
    );

    for (const entry& dEntry : layersDict)
    {
        if (dEntry.isDict())
        {
            const keyType& key = dEntry.keyword();
            const dictionary& layerDict = dEntry.dict();

            const labelHashSet patchIDs
            (
                boundaryMesh.patchSet(wordRes(one{}, wordRe(key)))
            );

            if (patchIDs.size() == 0)
            {
                IOWarningInFunction(layersDict)
                    << "Layer specification for " << key
                    << " does not match any patch." << endl
                    << "Valid patches are " << boundaryMesh.names() << endl;
            }
            else
            {
                for (const label patchi : patchIDs)
                {
                    numLayers_[patchi] =
                        layerDict.get<label>("nSurfaceLayers");

                    word spec;
                    if (layerDict.readIfPresent("thicknessModel", spec))
                    {
                        // If the thickness model is explicitly specified
                        // we want full specification of all parameters
                        layerModels_[patchi] = thicknessModelTypeNames_[spec];
                        readLayerParameters
                        (
                            false,                  // verbose
                            layerDict,
                            layerModels_[patchi],   // spec
                            firstLayerThickness_[patchi],
                            finalLayerThickness_[patchi],
                            thickness_[patchi],
                            expansionRatio_[patchi]
                        );
                        minThickness_[patchi] =
                            layerDict.get<scalar>("minThickness");
                    }
                    else
                    {
                        // Optional override of thickness parameters
                        switch (layerModels_[patchi])
                        {
                            case FIRST_AND_TOTAL:
                                layerDict.readIfPresent
                                (
                                    "firstLayerThickness",
                                    firstLayerThickness_[patchi]
                                );
                                layerDict.readIfPresent
                                (
                                    "thickness",
                                    thickness_[patchi]
                                );
                            break;

                            case FIRST_AND_EXPANSION:
                                layerDict.readIfPresent
                                (
                                    "firstLayerThickness",
                                    firstLayerThickness_[patchi]
                                );
                                layerDict.readIfPresent
                                (
                                    "expansionRatio",
                                    expansionRatio_[patchi]
                                );
                            break;

                            case FINAL_AND_TOTAL:
                                layerDict.readIfPresent
                                (
                                    "finalLayerThickness",
                                    finalLayerThickness_[patchi]
                                );
                                layerDict.readIfPresent
                                (
                                    "thickness",
                                    thickness_[patchi]
                                );
                            break;

                            case FINAL_AND_EXPANSION:
                                layerDict.readIfPresent
                                (
                                    "finalLayerThickness",
                                    finalLayerThickness_[patchi]
                                );
                                layerDict.readIfPresent
                                (
                                    "expansionRatio",
                                    expansionRatio_[patchi]
                                );
                            break;

                            case TOTAL_AND_EXPANSION:
                                layerDict.readIfPresent
                                (
                                    "thickness",
                                    thickness_[patchi]
                                );
                                layerDict.readIfPresent
                                (
                                    "expansionRatio",
                                    expansionRatio_[patchi]
                                );
                            break;

                            case FIRST_AND_RELATIVE_FINAL:
                                layerDict.readIfPresent
                                (
                                    "firstLayerThickness",
                                    firstLayerThickness_[patchi]
                                );
                                layerDict.readIfPresent
                                (
                                    "finalLayerThickness",
                                    finalLayerThickness_[patchi]
                                );
                            break;

                            default:
                                FatalIOErrorInFunction(dict)
                                    << "problem." << exit(FatalIOError);
                            break;
                        }

                        layerDict.readIfPresent
                        (
                            "minThickness",
                            minThickness_[patchi]
                        );
                    }

                    layerDict.readIfPresent
                    (
                        "relativeSizes",
                        relativeSizes_[patchi]
                    );
                }
            }
        }
    }


    forAll(numLayers_, patchi)
    {
        // Calculate the remaining parameters
        calculateLayerParameters
        (
            layerModels_[patchi],
            numLayers_[patchi],
            firstLayerThickness_[patchi],
            finalLayerThickness_[patchi],
            thickness_[patchi],
            expansionRatio_[patchi]
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::layerParameters::layerThickness
(
    const label nLayers,
    const scalar layerThickness,    // overall layer thickness
    const scalar expansionRatio,

    const label layerStart,         // start in nLayers
    const label layerSize           // size of slice of nLayers
)
{
    if (layerSize == 0 || nLayers == 0)
    {
        return 0.0;
    }
    else if (layerSize > nLayers || layerStart >= nLayers)
    {
        FatalErrorInFunction    << "Illegal input for slice of layer:"
            << " overall nLayers:" << nLayers
            << " slice nLayers:" << layerSize
            << " slice start:" << layerStart
            << exit(FatalError);
        return 0.0;
    }
    else if (mag(expansionRatio-1) < SMALL)
    {
        return layerThickness*layerSize/nLayers;
    }
    else
    {
        const scalar firstLayerThickness =
            finalLayerThicknessRatio(nLayers, expansionRatio)
          * layerThickness
          / pow(expansionRatio, nLayers-1);

        // Calculate thickness of single layer at layerStart
        const scalar startThickness =
            firstLayerThickness
           *pow(expansionRatio, layerStart);

        // See below for formula
        const scalar thickness =
            startThickness
           *(1.0 - pow(expansionRatio, layerSize))
           /(1.0 - expansionRatio);

        return thickness;
    }
}


Foam::scalar Foam::layerParameters::layerThickness
(
    const thicknessModelType layerSpec,
    const label nLayers,
    const scalar firstLayerThickness,
    const scalar finalLayerThickness,
    const scalar totalThickness,
    const scalar expansionRatio
)
{
    switch (layerSpec)
    {
        case FIRST_AND_TOTAL:
        case FINAL_AND_TOTAL:
        case TOTAL_AND_EXPANSION:
        {
            return totalThickness;
        }
        break;

        case FIRST_AND_EXPANSION:
        {
            if (mag(expansionRatio-1) < SMALL)
            {
                return firstLayerThickness * nLayers;
            }
            else
            {
                return firstLayerThickness
                   *(1.0 - pow(expansionRatio, nLayers))
                   /(1.0 - expansionRatio);
            }
        }
        break;

        case FINAL_AND_EXPANSION:
        {
            if (mag(expansionRatio-1) < SMALL)
            {
                return finalLayerThickness * nLayers;
            }
            else
            {
                scalar invExpansion = 1.0 / expansionRatio;
                return finalLayerThickness
                   *(1.0 - pow(invExpansion, nLayers))
                   /(1.0 - invExpansion);
            }
        }
        break;

        case FIRST_AND_RELATIVE_FINAL:
        {
            if (mag(expansionRatio-1) < SMALL)
            {
                return firstLayerThickness * nLayers;
            }
            else
            {
                scalar ratio = layerExpansionRatio
                (
                    layerSpec,
                    nLayers,
                    firstLayerThickness,
                    finalLayerThickness,
                    totalThickness,
                    expansionRatio
                );

                if (mag(ratio-1) < SMALL)
                {
                    return firstLayerThickness * nLayers;
                }
                else
                {
                    return firstLayerThickness *
                        (1.0 - pow(ratio, nLayers))
                      / (1.0 - ratio);
                }
            }
        }
        break;

        default:
        {
            FatalErrorInFunction
                << layerSpec << exit(FatalError);
            return -VGREAT;
        }
    }
}


Foam::scalar Foam::layerParameters::layerExpansionRatio
(
    const thicknessModelType layerSpec,
    const label nLayers,
    const scalar firstLayerThickness,
    const scalar finalLayerThickness,
    const scalar totalThickness,
    const scalar expansionRatio
)
{
    switch (layerSpec)
    {
        case FIRST_AND_EXPANSION:
        case FINAL_AND_EXPANSION:
        case TOTAL_AND_EXPANSION:
        {
            return expansionRatio;
        }
        break;

        case FIRST_AND_TOTAL:
        {
            if (firstLayerThickness < SMALL)
            {
                // Do what?
                return 1;
            }
            else
            {
                return layerExpansionRatio
                (
                    nLayers,
                    totalThickness/firstLayerThickness
                );
            }
        }
        break;

        case FINAL_AND_TOTAL:
        {
            if (finalLayerThickness < SMALL)
            {
                // Do what?
                return 1;
            }
            else
            {
                return
                    1.0
                  / layerExpansionRatio
                    (
                        nLayers,
                        totalThickness/finalLayerThickness
                    );
            }
        }
        break;

        case FIRST_AND_RELATIVE_FINAL:
        {
            if (firstLayerThickness < SMALL || nLayers <= 1)
            {
                return 1.0;
            }
            else
            {
                // Note: at this point the finalLayerThickness is already
                // absolute
                return pow
                (
                    finalLayerThickness/firstLayerThickness,
                    1.0/(nLayers-1)
                );
            }
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Illegal thickness specification" << exit(FatalError);
            return -VGREAT;
        }
    }
}


Foam::scalar Foam::layerParameters::firstLayerThickness
(
    const thicknessModelType layerSpec,
    const label nLayers,
    const scalar firstLayerThickness,
    const scalar finalLayerThickness,
    const scalar totalThickness,
    const scalar expansionRatio
)
{
    switch (layerSpec)
    {
        case FIRST_AND_EXPANSION:
        case FIRST_AND_TOTAL:
        case FIRST_AND_RELATIVE_FINAL:
        {
            return firstLayerThickness;
        }

        case FINAL_AND_EXPANSION:
        {
            if (expansionRatio < SMALL)
            {
                // Do what?
                return 0.0;
            }
            else
            {
                return finalLayerThickness*pow(1.0/expansionRatio, nLayers-1);
            }
        }
        break;

        case FINAL_AND_TOTAL:
        {
            scalar r = layerExpansionRatio
            (
                layerSpec,
                nLayers,
                firstLayerThickness,
                finalLayerThickness,
                totalThickness,
                expansionRatio
            );
            return finalLayerThickness/pow(r, nLayers-1);
        }
        break;

        case TOTAL_AND_EXPANSION:
        {
            scalar r = finalLayerThicknessRatio
            (
                nLayers,
                expansionRatio
            );
            scalar finalThickness = r*totalThickness;
            return finalThickness/pow(expansionRatio, nLayers-1);
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "Illegal thickness specification" << exit(FatalError);
            return -VGREAT;
        }
    }
}


Foam::scalar Foam::layerParameters::finalLayerThicknessRatio
(
    const label nLayers,
    const scalar expansionRatio
)
{
    if (nLayers > 0)
    {
        if (mag(expansionRatio-1) < SMALL)
        {
            return 1.0/nLayers;
        }
        else
        {
            return
                pow(expansionRatio, nLayers - 1)
               *(1.0 - expansionRatio)
               /(1.0 - pow(expansionRatio, nLayers));
        }
    }
    else
    {
        return 0.0;
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "gltfCoordSetWriter.H"
#include "coordSet.H"
#include "fileName.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "foamGltfScene.H"
#include "foamGltfSceneWriter.H"
#include "coordSetWriterMethods.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSetWriters
{
    defineTypeName(gltfWriter);
    addToRunTimeSelectionTable(coordSetWriter, gltfWriter, word);
    addToRunTimeSelectionTable(coordSetWriter, gltfWriter, wordDict);
}
}

const Foam::Enum<Foam::coordSetWriters::gltfWriter::fieldOption>
Foam::coordSetWriters::gltfWriter::fieldOptionNames_
({
    // No naming for NONE
    { fieldOption::UNIFORM, "uniform" },
    { fieldOption::FIELD, "field" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static tmp<vectorField> getBoundedColours
(
    const colourTable& colours,
    const Field<Type>& field,
    const scalar boundMin,
    const scalar boundMax
)
{
    const scalar boundDelta = (boundMax - boundMin + ROOTVSMALL);

    auto tresult = tmp<vectorField>::New(field.size());
    auto& result = tresult.ref();

    forAll(field, i)
    {
        const Type& val = field[i];

        const scalar f =
        (
            pTraits<Type>::nComponents == 1
          ? scalar(component(val, 0))
          : scalar(mag(val))
        );

        // 0-1 clipped by value()
        result[i] = colours.value((f - boundMin)/boundDelta);
    }

    return tresult;
}


template<class Type>
static vector getAnimationColour
(
    const dictionary& dict,
    const colourTable& colours,
    const Field<Type>& field
)
{
    scalar refValue(0);
    scalarMinMax valLimits;

    if (pTraits<Type>::nComponents == 1)
    {
        MinMax<Type> scanned(minMax(field));

        refValue = scalar(component(field[0], 0));
        valLimits.min() = scalar(component(scanned.min(), 0));
        valLimits.max() = scalar(component(scanned.max(), 0));
    }
    else
    {
        // Use mag() for multiple components
        refValue = mag(field[0]);
        valLimits = minMaxMag(field);
    }

    dict.readIfPresent("min", valLimits.min());
    dict.readIfPresent("max", valLimits.max());

    const scalar fraction =
    (
        (refValue - valLimits.min())
      / (valLimits.max() - valLimits.min() + ROOTVSMALL)
    );

    // 0-1 clipped by value()
    return colours.value(fraction);
}


} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::coordSetWriters::gltfWriter::getColourMap
(
    const dictionary& dict
) const
{
    word colourMap = colourTable::predefinedNames.names()[0];
    dict.readIfPresent("colourMap", colourMap);

    return colourMap;
}


const Foam::colourTable& Foam::coordSetWriters::gltfWriter::getColourTable
(
    const dictionary& dict
) const
{
    return colourTable::ref(getColourMap(dict));
}


Foam::scalarMinMax Foam::coordSetWriters::gltfWriter::getFieldLimits
(
    const word& fieldName
) const
{
    const dictionary fieldDict = fieldInfoDict_.subOrEmptyDict(fieldName);

    scalarMinMax limits(-GREAT, GREAT);

    fieldDict.readIfPresent("min", limits.min());
    fieldDict.readIfPresent("max", limits.max());

    return limits;
}


Foam::tmp<Foam::scalarField>
Foam::coordSetWriters::gltfWriter::getAlphaField
(
    const dictionary& dict
) const
{
    // Fallback value
    scalar alphaValue(1);

    const entry* eptr = dict.findEntry("alpha", keyType::LITERAL);

    if (!eptr)
    {
        // Not specified
    }
    else if (!eptr->stream().peek().isWord())
    {
        // Value specified

        ITstream& is = eptr->stream();
        is >> alphaValue;
        dict.checkITstream(is, "alpha");
    }
    else
    {
        // Enumeration

        const auto option = fieldOptionNames_.get("alpha", dict);

        switch (option)
        {
            case fieldOption::NONE:
            {
                break;
            }
            case fieldOption::UNIFORM:
            {
                dict.readEntry("alphaValue", alphaValue);
                break;
            }
            case fieldOption::FIELD:
            {
                WarningInFunction
                    << "Unsupported 'field' specification for alpha values"
                    << endl;
                break;
            }
        }
    }

    return tmp<scalarField>::New(1, alphaValue);
}


void Foam::coordSetWriters::gltfWriter::setupAnimationColour()
{
    const dictionary& dict = animationDict_;

    const entry* eptr = dict.findEntry("colour", keyType::LITERAL);

    if (!eptr || !eptr->isStream())
    {
        FatalIOErrorInFunction(dict)
            << "Missing 'colour' entry"
            << exit(FatalIOError);
    }
    else if (!eptr->stream().peek().isString())
    {
        // Value specified

        ITstream& is = eptr->stream();
        is >> animateColourValue_;
        dict.checkITstream(is, "colour");

        // Has uniform value
        animateColourOption_ = fieldOption::UNIFORM;
    }
    else
    {
        // Enumeration

        const auto option = fieldOptionNames_.get("colour", dict);

        switch (option)
        {
            case fieldOption::NONE:
            {
                FatalErrorInFunction
                    << "Cannot select 'none' for colour entry!" << nl
                    << "... possible programming error"
                    << exit(FatalError);
                break;
            }
            case fieldOption::UNIFORM:
            {
                dict.readEntry("colourValue", animateColourValue_);

                // Has uniform value
                animateColourOption_ = fieldOption::UNIFORM;
                break;
            }
            case fieldOption::FIELD:
            {
                // Needs named field...
                animateColourName_ = dict.get<word>("colourField");
                animateColourOption_ = fieldOption::FIELD;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSetWriters::gltfWriter::gltfWriter()
:
    coordSetWriter(),
    writer_(nullptr),
    animate_(false),
    colour_(false),
    animateColourOption_(fieldOption::NONE),
    animateColourName_(),
    animateColourValue_(Zero),
    fieldInfoDict_(),
    animationDict_()
{}


Foam::coordSetWriters::gltfWriter::gltfWriter(const dictionary& options)
:
    coordSetWriter(options),
    writer_(nullptr),
    animate_(options.getOrDefault("animate", false)),
    colour_(options.getOrDefault("colour", false)),
    animateColourOption_(fieldOption::NONE),
    animateColourName_(),
    animateColourValue_(Zero),
    fieldInfoDict_(options.subOrEmptyDict("fieldInfo")),
    animationDict_(options.subOrEmptyDict("animationInfo"))
{
    // fieldInfo
    // {
    //     U
    //     {
    //         colourMap       coolToWarm; // others...
    //         min             10;
    //         max             100;
    //         alpha           0.5;
    //     }
    // }
}


Foam::coordSetWriters::gltfWriter::gltfWriter
(
    const coordSet& coords,
    const fileName& outputPath,
    const dictionary& options
)
:
    gltfWriter(options)
{
    open(coords, outputPath);
}


Foam::coordSetWriters::gltfWriter::gltfWriter
(
    const UPtrList<coordSet>& tracks,
    const fileName& outputPath,
    const dictionary& options
)
:
    gltfWriter(options)
{
    open(tracks, outputPath);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coordSetWriters::gltfWriter::~gltfWriter()
{
    close();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::coordSetWriters::gltfWriter::path() const
{
    // 1) rootdir/<TIME>/setName.gltf
    // 2) rootdir/setName.gltf

    return getExpectedPath("gltf");
}


void Foam::coordSetWriters::gltfWriter::close(bool force)
{
    writer_.reset(nullptr);
    coordSetWriter::close(force);
}


void Foam::coordSetWriters::gltfWriter::beginTime(const Time& t)
{
    writer_.reset(nullptr);
    coordSetWriter::beginTime(t);
}


void Foam::coordSetWriters::gltfWriter::beginTime(const instant& inst)
{
    writer_.reset(nullptr);
    coordSetWriter::beginTime(inst);
}


void Foam::coordSetWriters::gltfWriter::endTime()
{
    writer_.reset(nullptr);
    coordSetWriter::endTime();
}


// * * * * * * * * * * * * * * * Implementation * * * * * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::coordSetWriters::gltfWriter::writeTemplate
(
    const word& fieldName,
    const UPtrList<const Field<Type>>& fieldPtrs
)
{
    if (coords_.size() != fieldPtrs.size())
    {
        FatalErrorInFunction
            << "Attempted to write field: " << fieldName
            << " (" << fieldPtrs.size() << " entries) for "
            << coords_.size() << " sets" << nl
            << exit(FatalError);
    }

    const auto& tracks = coords_;

    // const auto& times = trackTimes_;

    if (!writer_)
    {
        // Field:
        // 1) rootdir/<TIME>/setName.gltf
        // 2) rootdir/setName.gltf

        fileName outputFile = path();

        writer_.reset(new glTF::sceneWriter(outputFile));

        auto& scene = writer_->getScene();

        meshes_.resize(tracks.size());

        forAll(tracks, tracki)
        {
            word meshName("track:" + Foam::name(tracki));
            if (tracks.size() == 1)
            {
                meshName = "points";
            }

            meshes_[tracki] = scene.addMesh(tracks[tracki], meshName);
        }
    }


    auto& scene = writer_->getScene();

    forAll(tracks, tracki)
    {
        const label meshi = meshes_[tracki];
        const auto& field = fieldPtrs[tracki];

        scene.addFieldToMesh(field, fieldName, meshi);

        if (colour_)
        {
            const dictionary dict = fieldInfoDict_.subOrEmptyDict(fieldName);
            const auto& colours = getColourTable(dict);

            const auto talpha = getAlphaField(dict);
            const scalarField& alpha = talpha();

            const scalarMinMax valLimits = getFieldLimits(fieldName);

            scalarMinMax fldLimits;
            if (pTraits<Type>::nComponents == 1)
            {
                MinMax<Type> scanned(minMax(field));

                fldLimits.min() = scalar(component(scanned.min(), 0));
                fldLimits.max() = scalar(component(scanned.max(), 0));
            }
            else
            {
                // Use mag() for multiple components
                fldLimits = minMaxMag(field);
            }

            // Generated field colours
            vectorField fieldColour
            (
                getBoundedColours
                (
                    colours,
                    field,
                    max(fldLimits.min(), valLimits.min()),  // boundMin
                    min(fldLimits.max(), valLimits.max())   // boundMax
                )
            );

            scene.addColourToMesh
            (
                fieldColour,
                "Colour:" + fieldName,
                meshi,
                alpha
            );
        }
    }

    return writer_().path();
}


template<class Type>
Foam::fileName Foam::coordSetWriters::gltfWriter::writeTemplate_animate
(
    const word& fieldName,
    const UPtrList<const Field<Type>>& fieldPtrs
)
{
    if (coords_.size() != fieldPtrs.size())
    {
        FatalErrorInFunction
            << "Attempted to write field: " << fieldName
            << " (" << fieldPtrs.size() << " entries) for "
            << coords_.size() << " sets" << nl
            << exit(FatalError);
    }

    const auto& tracks = this->coords_;
    const auto& times = this->trackTimes_;

    if (!writer_)
    {
        // Field:
        // 1) rootdir/<TIME>/setName.gltf
        // 2) rootdir/setName.gltf

        fileName outputFile = path();

        writer_.reset(new glTF::sceneWriter(outputFile));

        auto& scene = writer_->getScene();

        meshes_.resize(tracks.size());

        const label animationi = scene.createAnimation("animation");

        forAll(tracks, tracki)
        {
            const auto& track = tracks[tracki];

            if (track.empty())
            {
                meshes_[tracki] = -1;
                continue;
            }

            // Seed starting position

            meshes_[tracki] =
                scene.addMesh
                (
                    vectorField(1, track[0]),
                    "track:" + Foam::name(tracki)
                );

            const label meshi = meshes_[tracki];

            // Time frames
            const label timeId =
                scene.addField(times[tracki], "time:" + Foam::name(tracki));

            // Translations
            const vectorField translation(track - track[0]);
            const label translationId =
                scene.addField(translation, "translation");

            scene.addToAnimation(animationi, timeId, translationId, meshi);
        }
    }


    auto& scene = writer_->getScene();

    // Seed starting field values

    forAll(tracks, tracki)
    {
        const auto& track = tracks[tracki];
        const label meshi = meshes_[tracki];
        const Field<Type>& field = fieldPtrs[tracki];

        if (track.empty() || meshi < 0)
        {
            continue;
        }

        // Seed starting field values
        scene.addFieldToMesh(Field<Type>(1, field[0]), fieldName, meshi);
    }


    // Note: colours cannot be animated... setting a fixed value.
    // However, we need to wait until the field is actually seen

    if (colour_)
    {
        if (animateColourOption_ == fieldOption::NONE)
        {
            // First time - scan for information
            setupAnimationColour();
        }

        switch (animateColourOption_)
        {
            case fieldOption::NONE:
            {
                // Should not occur
                break;
            }
            case fieldOption::UNIFORM:
            {
                // Colour value is known

                vectorField fieldColour(1, animateColourValue_);
                scalarField alphaChannel(1, 1.0);

                const auto talpha = getAlphaField(animationDict_);

                if (talpha && talpha().size())
                {
                    alphaChannel[0] = talpha()[0];
                }

                forAll(tracks, tracki)
                {
                    const auto& track = tracks[tracki];
                    const label meshi = meshes_[tracki];

                    if (track.empty() || meshi < 0)
                    {
                        continue;
                    }

                    scene.addColourToMesh
                    (
                        fieldColour,
                        "Colour:fixed",  // ... or "Colour:constant"
                        meshi,
                        alphaChannel
                    );
                }

                // Mark as done
                animateColourName_.clear();
                animateColourOption_ = fieldOption::FIELD;
                break;
            }
            case fieldOption::FIELD:
            {
                if
                (
                    !animateColourName_.empty()
                 && animateColourName_ == fieldName
                )
                {
                    // This is the desired colour field. Process now

                    const auto& colours = getColourTable(animationDict_);

                    vectorField fieldColour(1, Zero);
                    scalarField alphaChannel(1, 1.0);

                    const auto talpha = getAlphaField(animationDict_);

                    if (talpha && talpha().size())
                    {
                        alphaChannel[0] = talpha()[0];
                    }

                    forAll(tracks, tracki)
                    {
                        const auto& track = tracks[tracki];
                        const label meshi = meshes_[tracki];
                        const Field<Type>& field = fieldPtrs[tracki];

                        if (track.empty() || meshi < 0)
                        {
                            continue;
                        }

                        fieldColour[0] =
                            getAnimationColour(animationDict_, colours, field);

                        scene.addColourToMesh
                        (
                            fieldColour,
                            "Colour:fixed",  // ... or "Colour:constant"
                            meshi,
                            alphaChannel
                        );
                    }

                    // Mark colouring as done. Avoid retriggering
                    animateColourName_.clear();
                    animateColourOption_ = fieldOption::FIELD;
                }
                break;
            }
        }
    }

    return writer_().path();
}


template<class Type>
Foam::fileName Foam::coordSetWriters::gltfWriter::writeTemplate
(
    const word& fieldName,
    const Field<Type>& values
)
{
    checkOpen();
    if (coords_.empty())
    {
        return fileName::null;
    }

    UPtrList<const Field<Type>> fieldPtrs(repackageFields(values));
    return writeTemplate(fieldName, fieldPtrs);
}


template<class Type>
Foam::fileName Foam::coordSetWriters::gltfWriter::writeTemplate
(
    const word& fieldName,
    const List<Field<Type>>& fieldValues
)
{
    checkOpen();
    if (coords_.empty())
    {
        return fileName::null;
    }

    UPtrList<const Field<Type>> fieldPtrs(repackageFields(fieldValues));

    if (animate_ && trackTimes_.size() >= coords_.size())
    {
        return writeTemplate_animate(fieldName, fieldPtrs);
    }

    return writeTemplate(fieldName, fieldPtrs);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Field writing methods
defineCoordSetWriterWriteFields(Foam::coordSetWriters::gltfWriter);


// ************************************************************************* //

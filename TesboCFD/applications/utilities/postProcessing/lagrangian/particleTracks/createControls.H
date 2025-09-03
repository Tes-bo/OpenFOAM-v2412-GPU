// Read particleTrackProperties dictionary and extract values from it

const word dictName("particleTrackProperties");

#include "setConstantMeshDictionaryIO.H"

IOdictionary propsDict(dictIO);

const word cloudName(propsDict.get<word>("cloud"));

label sampleFrequency(propsDict.get<label>("sampleFrequency"));

label maxPositions(propsDict.get<label>("maxPositions"));

label maxTracks(propsDict.getOrDefault<label>("maxTracks", -1));

word setFormat(propsDict.getOrDefault<word>("setFormat", "vtk"));

// Optional - if empty, select all
wordRes acceptFields;
propsDict.readIfPresent("fields", acceptFields);

// Optional
wordRes excludeFields;
propsDict.readIfPresent("exclude", excludeFields);

const word UName(propsDict.getOrDefault<word>("U", "U"));

const dictionary formatOptions
(
    propsDict.subOrEmptyDict("formatOptions", keyType::LITERAL)
);


// ************************************************************************* //

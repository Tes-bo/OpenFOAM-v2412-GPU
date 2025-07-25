// Read particleTrackDict dictionary and extract values from it

const word dictName("particleTrackDict");

#include "setConstantMeshDictionaryIO.H"

IOdictionary propsDict(dictIO);

const word cloudName(propsDict.get<word>("cloud"));

// Mandatory - if empty, select none
wordRes acceptFields(propsDict.get<wordRes>("fields"));

// Optional
wordRes excludeFields;
propsDict.readIfPresent("exclude", excludeFields);


// ************************************************************************* //

// Force dlOpen of FOAM_DLOPEN_LIBS,
// which is a comma-separated list of quoted library names
// (eg, "finiteVolume","fvOptions","meshTools")
//
// Principally for Windows applications where library loading may otherwise
// be incomplete.

#ifdef FOAM_DLOPEN_LIBS

Foam::dlOpen({ FOAM_DLOPEN_LIBS }, false);  // Silent on errors

#endif

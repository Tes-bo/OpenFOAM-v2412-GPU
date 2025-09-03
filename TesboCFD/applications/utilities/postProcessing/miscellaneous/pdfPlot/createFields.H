    IOdictionary pdfDictionary
    (
        IOobject
        (
            "pdfDict",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const label nIntervals(pdfDictionary.get<label>("nIntervals"));

    const label nSamples(pdfDictionary.get<label>("nSamples"));

    const bool writeData(pdfDictionary.get<bool>("writeData"));


    const fileName pdfPath = runTime.path()/"pdf";
    mkDir(pdfPath);

    Random rndGen;

    autoPtr<distributionModel> p
    (
        distributionModel::New
        (
            pdfDictionary,
            rndGen
        )
    );

    const scalar xMin = p->minValue();
    const scalar xMax = p->maxValue();

    autoPtr<OFstream> filePtr;
    if (writeData)
    {
        fileName fName = pdfPath/(p->type() + ".data");
        Info<< "Writing " << p->type() << " data samples to file:" << nl
            << fName << nl << endl;

        filePtr.reset(new OFstream(fName));
    }

    scalarField samples(nIntervals, Zero);

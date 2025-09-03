// Maximum length for dynamicList
const label maxDynListLength
(
    viewFactorDict.getOrDefault<label>("maxDynListLength", 1000000000)
);

for (const int proci : Pstream::allProcs())
{
    std::vector<Point> start;
    start.reserve(coarseMesh.nFaces());

    std::vector<Point> end(start.size());
    end.reserve(start.size());

    DynamicList<label> startIndex(start.size());
    DynamicList<label> endIndex(start.size());

    DynamicList<label> startAgg(start.size());
    DynamicList<label> endAgg(start.size());

    const pointField& myFc = remoteCoarseCf[Pstream::myProcNo()];
    const vectorField& myArea = remoteCoarseSf[Pstream::myProcNo()];
    const labelField& myAgg = remoteCoarseAgg[Pstream::myProcNo()];

    const pointField& remoteArea = remoteCoarseSf[proci];
    const pointField& remoteFc = remoteCoarseCf[proci];
    const labelField& remoteAgg = remoteCoarseAgg[proci];

    label i = 0;
    label j = 0;
    do
    {
        for (; i < myFc.size(); i++)
        {
            const point& fc = myFc[i];
            const vector& fA = myArea[i];
            const label& fAgg = myAgg[i];

            for (; j < remoteFc.size(); j++)
            {
                if (proci != Pstream::myProcNo() || i != j)
                {
                    const point& remFc = remoteFc[j];
                    const vector& remA = remoteArea[j];
                    const label& remAgg = remoteAgg[j];

                    const vector d(remFc - fc);


                    const vector nd = d/mag(d);
                    const vector nfA = fA/mag(fA);
                    const vector nremA = remA/mag(remA);

                    if (((nd & nfA) < 0) && ((nd & nremA) > 0))
                    {
                        vector direction(d[0], d[1], d[2]);
                        vector s(fc[0], fc[1], fc[2]);
                        vector rayEnd(s + (1-intTol)*direction);
                        end.push_back(Point(rayEnd[0], rayEnd[1], rayEnd[2]));

                        s += vector(intTol*d[0], intTol*d[1], intTol*d[2]);

                        start.push_back(Point(s[0], s[1], s[2]));
                        startIndex.append(i);
                        if (useAgglomeration)
                        {
                            startAgg.append
                            (
                                globalNumbering.toGlobal(proci, fAgg)
                            );
                            endAgg.append
                            (
                                globalNumbering.toGlobal(proci, remAgg)
                            );
                        }

                        label globalI = globalNumbering.toGlobal(proci, j);
                        endIndex.append(globalI);

                        if (startIndex.size() > maxDynListLength)
                        {
                            FatalErrorInFunction
                                << "Dynamic list need from capacity."
                                << "Actual size maxDynListLength : "
                                <<  maxDynListLength
                                << abort(FatalError);
                        }
                    }
                }
            }

            if (j == remoteFc.size())
            {
                j = 0;
            }
        }

    } while (i < myFc.size());

    label totalRays(startIndex.size());
    reduce(totalRays, sumOp<scalar>());

    if (debug)
    {
        Pout<< "Number of rays :" << totalRays << endl;
    }

    for (unsigned long rayI = 0; rayI < start.size(); ++rayI)
    {
        Segment ray(start[rayI], end[rayI]);
        Segment_intersection intersects = tree.any_intersection(ray);

        if (!intersects)
        {
            rayStartFace.append(startIndex[rayI]);
            rayEndFace.append(endIndex[rayI]);
        }
    }
    start.clear();

    if (debug)
    {
        Pout << "hits : " <<  rayStartFace.size()<< endl;
    }

    startIndex.clear();
    end.clear();
    endIndex.clear();
    startAgg.clear();
    endAgg.clear();
}

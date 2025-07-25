// All rays expressed as start face (local) index end end face (global)
// Pre-size by assuming a certain percentage is visible.


if (Pstream::master())
{
    Info<< "\nShooting rays using distributedTriSurfaceMesh (no CGAL support)"
        << endl;
}



// Maximum length for dynamicList
const label maxDynListLength
(
    viewFactorDict.getOrDefault<label>("maxDynListLength", 1000000)
);

for (const int proci : Pstream::allProcs())
{
    // Shoot rays from me to proci. Note that even if processor has
    // 0 faces we still need to call findLine to keep calls synced.

    DynamicField<point> start(coarseMesh.nFaces());
    DynamicField<point> end(start.size());
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

                    if (((d & fA) < 0.) && ((d & remA) > 0))
                    {
                        start.append(fc + 0.001*d);
                        startIndex.append(i);
                        startAgg.append(globalNumbering.toGlobal(proci, fAgg));
                        end.append(fc + 0.999*d);
                        label globalI = globalNumbering.toGlobal(proci, j);
                        endIndex.append(globalI);
                        endAgg.append(globalNumbering.toGlobal(proci, remAgg));
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

    } while (returnReduceOr(i < myFc.size()));

    List<pointIndexHit> hitInfo(startIndex.size());
    surfacesMesh.findLine(start, end, hitInfo);

    // Return hit global agglo index
    labelList aggHitIndex;
    surfacesMesh.getField(hitInfo, aggHitIndex);

    DynamicList<label> dRayIs;

    // Collect the rays which have no obstacle in between rayStartFace
    // and rayEndFace. If the ray hit itself, it gets stored in dRayIs
    forAll(hitInfo, rayI)
    {
        if (!hitInfo[rayI].hit())
        {
            rayStartFace.append(startIndex[rayI]);
            rayEndFace.append(endIndex[rayI]);
        }
        else if (aggHitIndex[rayI] == startAgg[rayI])
        {
            dRayIs.append(rayI);
        }
    }

    start.clear();


    // Continue rays which hit themself. If they hit the target
    // agglomeration are added to rayStartFace and rayEndFace

    bool firstLoop = true;
    DynamicField<point> startHitItself;
    DynamicField<point> endHitItself;
    label iter = 0;

    do
    {
        labelField rayIs;
        rayIs.transfer(dRayIs);
        dRayIs.clear();
        forAll(rayIs, rayI)
        {
            const label rayID = rayIs[rayI];
            label hitIndex = -1;

            if (firstLoop)
            {
                hitIndex = rayIs[rayI];
            }
            else
            {
                hitIndex = rayI;
            }

            if (hitInfo[hitIndex].hit())
            {
                if (aggHitIndex[hitIndex] == startAgg[rayID])
                {
                    const vector& endP = end[rayID];
                    const vector& startP = hitInfo[hitIndex].point();
                    const vector& d = endP - startP;

                    startHitItself.append(startP + 0.01*d);
                    endHitItself.append(startP + 1.01*d);

                    dRayIs.append(rayID);
                }
                else if (aggHitIndex[hitIndex] == endAgg[rayID])
                {
                    rayStartFace.append(startIndex[rayID]);
                    rayEndFace.append(endIndex[rayID]);
                }

            }
        }

        hitInfo.clear();
        hitInfo.resize(dRayIs.size());

        surfacesMesh.findLine(startHitItself, endHitItself, hitInfo);

        surfacesMesh.getField(hitInfo, aggHitIndex);


        endHitItself.clear();
        startHitItself.clear();
        firstLoop = false;
        iter ++;

    } while (returnReduceOr(hitInfo.size()) && iter < 10);

    startIndex.clear();
    end.clear();
    endIndex.clear();
    startAgg.clear();
    endAgg.clear();
}

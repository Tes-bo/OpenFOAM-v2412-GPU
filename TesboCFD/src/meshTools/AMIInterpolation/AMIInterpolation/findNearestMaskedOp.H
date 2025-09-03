#include "indexedOctree.H"
#include "treeDataPrimitivePatch.H"

namespace Foam
{

template<class PatchType>
class findNearestMaskedOp
{
    const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree_;
    const labelUList& excludeIndices_;

public:

    findNearestMaskedOp
    (
        const indexedOctree<treeDataPrimitivePatch<PatchType>>& tree,
        const labelUList& excludeIndices
    )
    :
        tree_(tree),
        excludeIndices_(excludeIndices)
    {}

    void operator()
    (
        const labelUList& indices,
        const point& sample,

        scalar& nearestDistSqr,
        label& minIndex,
        point& nearestPoint
    ) const
    {
        const treeDataPrimitivePatch<PatchType>& shape = tree_.shapes();
        const PatchType& patch = shape.patch();

        const pointField& points = patch.points();

        forAll(indices, i)
        {
            const label index = indices[i];

            if (!excludeIndices_.found(index))
            {
                const auto& f = patch[index];

                pointHit nearHit = f.nearestPoint(sample, points);
                scalar distSqr = sqr(nearHit.distance());

                if (distSqr < nearestDistSqr)
                {
                    nearestDistSqr = distSqr;
                    minIndex = index;
                    nearestPoint = nearHit.point();
                }
            }
        }
    }
};

} // End namespace Foam

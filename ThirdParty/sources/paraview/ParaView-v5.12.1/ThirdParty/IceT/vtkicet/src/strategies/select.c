/* -*- c -*- *******************************************************/
/*
 * Copyright (C) 2010 Sandia Corporation
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * This source code is released under the New BSD License.
 */

#include <IceTDevStrategySelect.h>

#include <IceT.h>
#include <IceTDevDiagnostics.h>
#include <IceTDevImage.h>
#include <IceTDevState.h>

/* Declaration of strategy compose functions. */
extern IceTImage icetDirectCompose(void);
extern IceTImage icetSequentialCompose(void);
extern IceTImage icetSplitCompose(void);
extern IceTImage icetReduceCompose(void);
extern IceTImage icetVtreeCompose(void);

/* Declaration of single image strategy compose functions. */
extern void icetAutomaticCompose(const IceTInt *compose_group,
                                 IceTInt group_size,
                                 IceTInt image_dest,
                                 IceTSparseImage input_image,
                                 IceTSparseImage *result_image,
                                 IceTSizeType *piece_offset);
extern void icetBswapCompose(const IceTInt *compose_group,
                             IceTInt group_size,
                             IceTInt image_dest,
                             IceTSparseImage input_image,
                             IceTSparseImage *result_image,
                             IceTSizeType *piece_offset);
void icetBswapFoldingCompose(const IceTInt *compose_group,
                             IceTInt group_size,
                             IceTInt image_dest,
                             IceTSparseImage input_image,
                             IceTSparseImage *result_image,
                             IceTSizeType *piece_offset);
extern void icetTreeCompose(const IceTInt *compose_group,
                            IceTInt group_size,
                            IceTInt image_dest,
                            IceTSparseImage input_image,
                            IceTSparseImage *result_image,
                            IceTSizeType *piece_offset);
extern void icetRadixkCompose(const IceTInt *compose_group,
                              IceTInt group_size,
                              IceTInt image_dest,
                              IceTSparseImage input_image,
                              IceTSparseImage *result_image,
                              IceTSizeType *piece_offset);
extern void icetRadixkrCompose(const IceTInt *compose_group,
                               IceTInt group_size,
                               IceTInt image_dest,
                               IceTSparseImage input_image,
                               IceTSparseImage *result_image,
                               IceTSizeType *piece_offset);

/*==================================================================*/

IceTBoolean icetStrategyValid(IceTEnum strategy)
{
    switch (strategy) {
      case ICET_STRATEGY_DIRECT:
      case ICET_STRATEGY_SEQUENTIAL:
      case ICET_STRATEGY_SPLIT:
      case ICET_STRATEGY_REDUCE:
      case ICET_STRATEGY_VTREE:
          return ICET_TRUE;
      default:
          return ICET_FALSE;
    }
}

const char *icetStrategyNameFromEnum(IceTEnum strategy)
{
    switch (strategy) {
      case ICET_STRATEGY_DIRECT:        return "Direct";
      case ICET_STRATEGY_SEQUENTIAL:    return "Sequential";
      case ICET_STRATEGY_SPLIT:         return "Split";
      case ICET_STRATEGY_REDUCE:        return "Reduce";
      case ICET_STRATEGY_VTREE:         return "Virtual Tree";
      case ICET_STRATEGY_UNDEFINED:
          icetRaiseError(ICET_INVALID_ENUM,
                         "Strategy not defined. "
                         "Use icetStrategy to set the strategy.");
          return "<Not Set>";
      default:
          icetRaiseError(ICET_INVALID_ENUM,
                         "Invalid strategy %d.", strategy);
          return "<Invalid>";
    }
}

IceTBoolean icetStrategySupportsOrdering(IceTEnum strategy)
{
    switch (strategy) {
      case ICET_STRATEGY_DIRECT:        return ICET_TRUE;
      case ICET_STRATEGY_SEQUENTIAL:    return ICET_TRUE;
      case ICET_STRATEGY_SPLIT:         return ICET_FALSE;
      case ICET_STRATEGY_REDUCE:        return ICET_TRUE;
      case ICET_STRATEGY_VTREE:         return ICET_FALSE;
      case ICET_STRATEGY_UNDEFINED:
          icetRaiseError(ICET_INVALID_ENUM,
                         "Strategy not defined. "
                         "Use icetStrategy to set the strategy.");
          return ICET_FALSE;
      default:
          icetRaiseError(ICET_INVALID_ENUM, "Invalid strategy %d.", strategy);
          return ICET_FALSE;
    }
}

IceTImage icetInvokeStrategy(IceTEnum strategy)
{
    icetRaiseDebug("Invoking strategy %s",
                   icetStrategyNameFromEnum(strategy));

    switch (strategy) {
      case ICET_STRATEGY_DIRECT:        return icetDirectCompose();
      case ICET_STRATEGY_SEQUENTIAL:    return icetSequentialCompose();
      case ICET_STRATEGY_SPLIT:         return icetSplitCompose();
      case ICET_STRATEGY_REDUCE:        return icetReduceCompose();
      case ICET_STRATEGY_VTREE:         return icetVtreeCompose();
      case ICET_STRATEGY_UNDEFINED:
          icetRaiseError(ICET_INVALID_ENUM,
                         "Strategy not defined. "
                         "Use icetStrategy to set the strategy.");
          return icetImageNull();
      default:
          icetRaiseError(ICET_INVALID_ENUM, "Invalid strategy %d.", strategy);
          return icetImageNull();
    }
}

/*==================================================================*/

IceTBoolean icetSingleImageStrategyValid(IceTEnum strategy)
{
    switch(strategy) {
      case ICET_SINGLE_IMAGE_STRATEGY_AUTOMATIC:
      case ICET_SINGLE_IMAGE_STRATEGY_BSWAP:
      case ICET_SINGLE_IMAGE_STRATEGY_TREE:
      case ICET_SINGLE_IMAGE_STRATEGY_RADIXK:
      case ICET_SINGLE_IMAGE_STRATEGY_RADIXKR:
      case ICET_SINGLE_IMAGE_STRATEGY_BSWAP_FOLDING:
          return ICET_TRUE;
      default:
          return ICET_FALSE;
    }
}

const char *icetSingleImageStrategyNameFromEnum(IceTEnum strategy)
{
    switch(strategy) {
      case ICET_SINGLE_IMAGE_STRATEGY_AUTOMATIC:        return "Automatic";
      case ICET_SINGLE_IMAGE_STRATEGY_BSWAP:            return "Binary Swap";
      case ICET_SINGLE_IMAGE_STRATEGY_TREE:             return "Binary Tree";
      case ICET_SINGLE_IMAGE_STRATEGY_RADIXK:           return "Radix-k";
      case ICET_SINGLE_IMAGE_STRATEGY_RADIXKR:          return "Radix-kr";
      case ICET_SINGLE_IMAGE_STRATEGY_BSWAP_FOLDING:    return "Folded Binary Swap";
      default:
          icetRaiseError(ICET_INVALID_ENUM,
                         "Invalid single image strategy %d.", strategy);
          return "<Invalid>";
    }
}

void icetInvokeSingleImageStrategy(IceTEnum strategy,
                                   const IceTInt *compose_group,
                                   IceTInt group_size,
                                   IceTInt image_dest,
                                   IceTSparseImage input_image,
                                   IceTSparseImage *result_image,
                                   IceTSizeType *piece_offset)
{
    icetRaiseDebug("Invoking single image strategy %s",
                   icetSingleImageStrategyNameFromEnum(strategy));

    switch(strategy) {
      case ICET_SINGLE_IMAGE_STRATEGY_AUTOMATIC:
          icetAutomaticCompose(compose_group,
                               group_size,
                               image_dest,
                               input_image,
                               result_image,
                               piece_offset);
          break;
      case ICET_SINGLE_IMAGE_STRATEGY_BSWAP:
          icetBswapCompose(compose_group,
                           group_size,
                           image_dest,
                           input_image,
                           result_image,
                           piece_offset);
          break;
      case ICET_SINGLE_IMAGE_STRATEGY_TREE:
          icetTreeCompose(compose_group,
                          group_size,
                          image_dest,
                          input_image,
                          result_image,
                          piece_offset);
          break;
      case ICET_SINGLE_IMAGE_STRATEGY_RADIXK:
          icetRadixkCompose(compose_group,
                            group_size,
                            image_dest,
                            input_image,
                            result_image,
                            piece_offset);
          break;
      case ICET_SINGLE_IMAGE_STRATEGY_RADIXKR:
          icetRadixkrCompose(compose_group,
                             group_size,
                             image_dest,
                             input_image,
                             result_image,
                             piece_offset);
          break;
    case ICET_SINGLE_IMAGE_STRATEGY_BSWAP_FOLDING:
        icetBswapFoldingCompose(compose_group,
                                group_size,
                                image_dest,
                                input_image,
                                result_image,
                                piece_offset);
        break;
      default:
          icetRaiseError(ICET_INVALID_ENUM,
                         "Invalid single image strategy %d.", strategy);
          break;
    }

    icetStateCheckMemory();
}

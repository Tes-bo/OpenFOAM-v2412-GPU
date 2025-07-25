/*
 * Copyright (c) 2013-2018 Intel, Inc. All rights reserved
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "mtl_ofi.h"

__opal_attribute_always_inline__ static inline int
ompi_mtl_ofi_irecv_false(struct mca_mtl_base_module_t *mtl,
               struct ompi_communicator_t *comm,
               int src,
               int tag,
               struct opal_convertor_t *convertor,
               mca_mtl_request_t *mtl_request)
{
    const bool OFI_CQ_DATA = false;

    return ompi_mtl_ofi_irecv_generic(mtl, comm, src, tag,
                                    convertor, mtl_request,
                                    OFI_CQ_DATA);
}

__opal_attribute_always_inline__ static inline int
ompi_mtl_ofi_irecv_true(struct mca_mtl_base_module_t *mtl,
               struct ompi_communicator_t *comm,
               int src,
               int tag,
               struct opal_convertor_t *convertor,
               mca_mtl_request_t *mtl_request)
{
    const bool OFI_CQ_DATA = true;

    return ompi_mtl_ofi_irecv_generic(mtl, comm, src, tag,
                                    convertor, mtl_request,
                                    OFI_CQ_DATA);
}

void ompi_mtl_ofi_irecv_symtable_init(struct ompi_mtl_ofi_symtable *sym_table)
{

    sym_table->ompi_mtl_ofi_irecv[false]
        = ompi_mtl_ofi_irecv_false;


    sym_table->ompi_mtl_ofi_irecv[true]
        = ompi_mtl_ofi_irecv_true;

}


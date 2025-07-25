/* -*- c -*- *****************************************************************
** Copyright (C) 2003 Sandia Corporation
** Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
** the U.S. Government retains certain rights in this software.
**
** This source code is released under the New BSD License.
**
** This test checks to make sure the size of compressed images never
** exceeds the advertised maximum buffer size.
*****************************************************************************/

#include "test_codes.h"
#include "test_util.h"

#include <IceTDevImage.h>
#include <IceTDevPorting.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

static void dump_image(const char *name,
                       IceTImage image)
{
    IceTInt rank;
    char filename[FILENAME_MAX];

    icetGetIntegerv(ICET_RANK, &rank);
    icetSnprintf(filename, FILENAME_MAX, "%s%03d.ppm", name, rank);
    write_ppm(filename,
              icetImageGetColorcub(image),
              (int)icetImageGetWidth(image),
              (int)icetImageGetHeight(image));
}

#define DIFF(x, y)      ((x) < (y) ? (y) - (x) : (x) - (y))

static int compare_ubyte_buffers(const IceTUByte *buffer1,
                                 const IceTUByte *buffer2,
                                 IceTSizeType buffer_size)
{
    IceTSizeType i;

    printstat("Checking ubyte buffers.\n");
    for (i = 0; i < buffer_size; i++) {
        if (DIFF(buffer1[i], buffer2[i]) > 5) {
            printstat("Bad pixel at %d!!!!!!!!!!\n", i);
            printstat("Got %d. Expected %d\n", buffer1[i], buffer2[i]);
            return ICET_FALSE;
        }
    }
    return ICET_TRUE;
}

static int compare_float_buffers(const IceTFloat *buffer1,
                                 const IceTFloat *buffer2,
                                 IceTSizeType buffer_size)
{
    IceTSizeType i;

    printstat("Checking float buffers.\n");
    for (i = 0; i < buffer_size; i++) {
        if (DIFF(buffer1[i], buffer2[i]) > 0.02f) {
            printstat("Bad pixel at %d!!!!!!!!!!\n", i);
            printstat("Got %f. Expected %f\n", buffer1[i], buffer2[i]);
            return ICET_FALSE;
        }
    }
    return ICET_TRUE;
}

static int compare_color_buffers(IceTImage refimage, IceTImage testimage) {
    int result;
    IceTSizeType num_pixels;

    if (   (icetImageGetWidth(refimage) != icetImageGetWidth(testimage))
        || (icetImageGetHeight(refimage) != icetImageGetHeight(testimage))) {
        printrank("Image dimensions not what is expected!!!!!!\n");
        printrank("Expected %dx%d, received %dx%d\n",
                  (int)icetImageGetWidth(refimage),
                  (int)icetImageGetHeight(refimage),
                  (int)icetImageGetWidth(testimage),
                  (int)icetImageGetHeight(testimage));
        return ICET_FALSE;
    }
    num_pixels = icetImageGetNumPixels(refimage);
    if (   icetImageGetColorFormat(refimage)
        != icetImageGetColorFormat(testimage)) {
        printrank("Color buffer formats do not match!!!!!!!\n");
        return ICET_FALSE;
    }

    switch (icetImageGetColorFormat(refimage)) {
    case ICET_IMAGE_COLOR_RGBA_UBYTE:
        result = compare_ubyte_buffers(icetImageGetColorub(refimage),
                                       icetImageGetColorub(testimage),
                                       num_pixels*4);
        break;
    case ICET_IMAGE_COLOR_RGBA_FLOAT:
        result = compare_float_buffers(icetImageGetColorcf(refimage),
                                       icetImageGetColorcf(testimage),
                                       num_pixels*4);
        break;
    case ICET_IMAGE_COLOR_RGB_FLOAT:
        result = compare_float_buffers(icetImageGetColorcf(refimage),
                                       icetImageGetColorcf(testimage),
                                       num_pixels*3);
        break;
    default:
        printrank("Unknown pixel format?!!!!!\n");
        return ICET_FALSE;
    }

    if (!result) {
        dump_image("bad", testimage);
        dump_image("ref", refimage);
    }

    return result;
}

static void check_image(IceTImage reference,
                        IceTImage image,
                        int *result)
{
    if (icetImageGetColorFormat(image) != ICET_IMAGE_COLOR_NONE) {
        if (!compare_color_buffers(reference, image)) {
            *result = TEST_FAILED;
        }
    }
}

static void InitImage(IceTImage image)
{
    /* Create a consistent test image in the given format. */
    IceTEnum format;
    IceTSizeType width;
    IceTSizeType height;
    IceTSizeType scale;

    width = icetImageGetWidth(image);
    height = icetImageGetHeight(image);
    scale = width + height;

    format = icetImageGetColorFormat(image);
    if (format == ICET_IMAGE_COLOR_RGBA_UBYTE) {
        IceTUByte *buffer = icetImageGetColorub(image);
        IceTSizeType x, y;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                buffer[0] = (IceTUByte)((255*x)/scale);
                buffer[1] = (IceTUByte)((255*y)/scale);
                buffer[2] = (IceTUByte)((255*(x+y))/scale);
                buffer[3] = (IceTUByte)255;
                buffer += 4;
            }
        }
    } else if (format == ICET_IMAGE_COLOR_RGBA_FLOAT) {
        IceTFloat *buffer = icetImageGetColorf(image);
        IceTSizeType x, y;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                buffer[0] = ((IceTFloat)x)/scale;
                buffer[1] = ((IceTFloat)y)/scale;
                buffer[2] = ((IceTFloat)(x+y))/scale;
                buffer[3] = (IceTFloat)1;
                buffer += 4;
            }
        }
    } else if (format == ICET_IMAGE_COLOR_RGB_FLOAT) {
        IceTFloat *buffer = icetImageGetColorf(image);
        IceTSizeType x, y;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                buffer[0] = ((IceTFloat)x)/scale;
                buffer[1] = ((IceTFloat)y)/scale;
                buffer[2] = ((IceTFloat)(x+y))/scale;
                buffer += 3;
            }
        }
    } else {
        printrank("*** Unknown color format? ***\n");
    }
}

static void ConvertImage(IceTImage source_image, IceTImage dest_image)
{
    switch (icetImageGetColorFormat(dest_image)) {
    case ICET_IMAGE_COLOR_RGBA_UBYTE:
        icetImageCopyColorub(source_image,
                             icetImageGetColorub(dest_image),
                             icetImageGetColorFormat(dest_image));
        break;
    case ICET_IMAGE_COLOR_RGBA_FLOAT:
    case ICET_IMAGE_COLOR_RGB_FLOAT:
        icetImageCopyColorf(source_image,
                            icetImageGetColorf(dest_image),
                            icetImageGetColorFormat(dest_image));
        break;
    case ICET_IMAGE_COLOR_NONE:
        /* Nothing to do. */
        break;
    default:
        printrank("Encountered unknown color format!\n");
    }

    switch (icetImageGetDepthFormat(dest_image)) {
    case ICET_IMAGE_DEPTH_FLOAT:
        icetImageCopyDepthf(source_image,
                            icetImageGetDepthf(dest_image),
                            icetImageGetDepthFormat(dest_image));
        break;
    case ICET_IMAGE_DEPTH_NONE:
        /* Nothing to do. */
        break;
    default:
        printrank("Encountered unknown color format!\n");
    }
}

static int DoImageConvertTest(IceTEnum src_color_format,
                              IceTEnum dest_color_format)
{
    IceTImage src_image;
    IceTVoid *src_buffer;
    IceTImage dest_image;
    IceTVoid *dest_buffer;
    IceTImage ref_image;
    IceTVoid *ref_buffer;
    IceTSizeType src_buffer_size;
    IceTSizeType dest_buffer_size;
    IceTEnum src_depth_format = ICET_IMAGE_DEPTH_NONE;
    IceTEnum dest_depth_format = ICET_IMAGE_DEPTH_NONE;
    int result;

    printstat("Using source color format of 0x%x\n", (int)src_color_format);
    printstat("Using dest color format of 0x%x\n", (int)dest_color_format);

    printstat("Allocating memory for %dx%d pixel image.\n",
              (int)SCREEN_WIDTH, (int)SCREEN_HEIGHT);

    icetSetColorFormat(src_color_format);
    icetSetDepthFormat(src_depth_format);
    src_buffer_size = icetImageBufferSize(SCREEN_WIDTH, SCREEN_HEIGHT);

    src_buffer = malloc(src_buffer_size);
    src_image = icetImageAssignBuffer(src_buffer, SCREEN_WIDTH, SCREEN_HEIGHT);

    icetSetColorFormat(dest_color_format);
    icetSetDepthFormat(dest_depth_format);
    dest_buffer_size = icetImageBufferSize(SCREEN_WIDTH, SCREEN_HEIGHT);

    dest_buffer = malloc(dest_buffer_size);
    dest_image = icetImageAssignBuffer(dest_buffer,SCREEN_WIDTH,SCREEN_HEIGHT);

    ref_buffer = malloc(dest_buffer_size);
    ref_image = icetImageAssignBuffer(ref_buffer, SCREEN_WIDTH, SCREEN_HEIGHT);

    printstat("Initializing images\n");
    InitImage(src_image);
    InitImage(ref_image);

    printstat("Doing conversion\n");
    ConvertImage(src_image, dest_image);

    printstat("Comparing conversion result\n");
    result = TEST_PASSED;
    if (!compare_color_buffers(ref_image, dest_image)) {
        result = TEST_FAILED;
    }

    printstat("Cleaning up.\n");
    free(src_buffer);
    free(dest_buffer);
    free(ref_buffer);

    return result;
}

#define NUM_COLOR_FORMATS 3
static IceTEnum COLOR_FORMATS[NUM_COLOR_FORMATS] = {
    ICET_IMAGE_COLOR_RGBA_UBYTE,
    ICET_IMAGE_COLOR_RGBA_FLOAT,
    ICET_IMAGE_COLOR_RGB_FLOAT
};

static int ImageConvertRun()
{
    int src_color_format_index;
    int dest_color_format_index;

    for (src_color_format_index = 0;
         src_color_format_index < NUM_COLOR_FORMATS;
         src_color_format_index++) {
        for (dest_color_format_index = 0;
             dest_color_format_index < NUM_COLOR_FORMATS;
             dest_color_format_index++) {
            int result =
                DoImageConvertTest(COLOR_FORMATS[src_color_format_index],
                                   COLOR_FORMATS[dest_color_format_index]);
            if (result != TEST_PASSED) {
                return result;
            }
        }
    }

    return TEST_PASSED;
}

int ImageConvert(int argc, char *argv[])
{
    /* To remove warning */
    (void)argc;
    (void)argv;

    return run_test(ImageConvertRun);
}

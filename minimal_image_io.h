#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <iostream>
#include <iostream>
#include <cmath>

// Reads an RGB image and returns a triple-pointer image_data[x][y][c],
// with x in [0,width-1], y in [0,height-1], c in [0,2], and pixel values in [0,255].
// Allocates memory that must be freed by the caller (both block[] and image_data[][][]).
bool read_image_as_3d_uchar(
    const char* filename,
    int &width, int &height,
    unsigned char ****image_data // (out) triple pointer
)
{
    int channels_in_file;
    // Force load as 3 channels (RGB):
    unsigned char* data = stbi_load(filename, &width, &height, &channels_in_file, 3);
    if (!data) {
        std::cerr << "Failed to load image: " << filename << std::endl;
        return false;
    }

    // We'll copy the image data from 'data' into our own allocated block,
    // so that we can safely free 'data' (stbi's allocation).
    int size = width * height * 3;
    unsigned char* block = new unsigned char[size];

    for (int i = 0; i < size; i++) {
        block[i] = data[i];
    }
    stbi_image_free(data); // We no longer need the stb buffer.

    // Allocate the top-level x dimension:
    *image_data = new unsigned char**[width];

    // Allocate an array for each x that can hold height pointers:
    for (int x = 0; x < width; x++) {
        (*image_data)[x] = new unsigned char*[height];
    }

    // Point each image_data[x][y] to the correct location in 'block'.
    // The pixel at (x,y) is at index = (y*width + x)*3 in block.
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            (*image_data)[x][y] = &block[(y*width + x)*3];
        }
    }

    return true;
}

bool write_image_as_3d_uchar(
    const char* filename,
    int width, int height,
    unsigned char ***image_data
)
{
    // We'll pack the triple-pointer data into a single contiguous array
    // to pass to stbi_write_png.
    unsigned char* out_data = new unsigned char[width * height * 3];

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            int idx = (y*width + x)*3;
            out_data[idx + 0] = image_data[x][y][0];
            out_data[idx + 1] = image_data[x][y][1];
            out_data[idx + 2] = image_data[x][y][2];
        }
    }

    // Write PNG (3 channels, stride is width*3 bytes per row)
    int stride_in_bytes = width * 3;
    bool success = (stbi_write_png(filename, width, height, 3, out_data, stride_in_bytes) != 0);
    if (!success) {
        std::cerr << "Failed to write image: " << filename << std::endl;
    }

    delete[] out_data;
    return success;
}

// Writes an RGB image from a triple-pointer image_data[x][y][c] in [0,1] to a PNG file.
// image_data must be width by height by 3.
bool write_image_as_3d_double(
    const char* filename,
    int width, int height,
    double ***image_data
)
{
    // Convert back to unsigned char
    unsigned char* out_data = new unsigned char[width * height * 3];

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            int idx = (y*width + x)*3;
            for (int c = 0; c < 3; c++) {
                double val = image_data[x][y][c] * 255.0;
                if (val < 0) val = 0;
                if (val > 255) val = 255;
                out_data[idx + c] = (unsigned char)(val);
            }
        }
    }

    // Write PNG
    int stride_in_bytes = width * 3;
    if (!stbi_write_png(filename, width, height, 3, out_data, stride_in_bytes)) {
        std::cerr << "Failed to write image: " << filename << std::endl;
        delete[] out_data;
        return false;
    }

    delete[] out_data;
    return true;
}
/*
==========================================================================
File:        image_processing.c 
Authors:     Costanza Maria Improta
==========================================================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jpeglib.h>


/** Read the JPEG image at `filename` as an array of bytes.
  Data is returned through the out pointers, while the return
  value indicates success or failure.
  NOTE: 1) if image is RGB, then the bytes are concatenated in R-G-B order
        2) `image` should be freed by the user
 */
static inline int
read_JPEG_file(char *filename,
               int *width, int *height, int *channels, unsigned char *(image[]))
{
  FILE *infile;
  if ((infile = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    return 0;
  }

  struct jpeg_error_mgr jerr;
  struct jpeg_decompress_struct cinfo;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, infile);
  (void) jpeg_read_header(&cinfo, TRUE);
  (void) jpeg_start_decompress(&cinfo);

  *width = cinfo.output_width, *height = cinfo.output_height;
  *channels = cinfo.num_components;
  // printf("width=%d height=%d c=%d\n", *width, *height, *channels);
  *image = malloc(*width * *height * *channels * sizeof(*image));
  JSAMPROW rowptr[1];
  int row_stride = *width * *channels;

  while (cinfo.output_scanline < cinfo.output_height) {
    rowptr[0] = *image + row_stride * cinfo.output_scanline;
    jpeg_read_scanlines(&cinfo, rowptr, 1);
  }
  jpeg_finish_decompress(&cinfo);

  jpeg_destroy_decompress(&cinfo);
  fclose(infile);
  return 1;
}


/** Writes the image in the specified file.
  NOTE: works with Grayscale or RGB modes only (based on number of channels)
 */
static inline void
write_JPEG_file(char *filename, int width, int height, int channels,
                unsigned char image[], int quality)
{
  FILE *outfile;
  if ((outfile = fopen(filename, "wb")) == NULL) {
    fprintf(stderr, "can't open %s\n", filename);
    exit(1);
  }

  struct jpeg_error_mgr jerr;
  struct jpeg_compress_struct cinfo;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo,outfile);

  cinfo.image_width = width;
  cinfo.image_height = height;
  cinfo.input_components = channels;
  cinfo.in_color_space = channels == 1 ? JCS_GRAYSCALE : JCS_RGB;
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, quality, TRUE);

  jpeg_start_compress(&cinfo, TRUE);
  JSAMPROW rowptr[1];
  int row_stride = width * channels;
  while (cinfo.next_scanline < cinfo.image_height) {
    rowptr[0] = & image[cinfo.next_scanline * row_stride];
    jpeg_write_scanlines(&cinfo, rowptr, 1);
  }
  jpeg_finish_compress(&cinfo);

  fclose(outfile);
  jpeg_destroy_compress(&cinfo);
}

/*
  The following method gets the minimum pixel value of the inputted image
*/
int getMin(unsigned char image[], int width, int height){
  int min = image[0];
  int index;
  for(index = 1; index < height*width; index++){
    if(min > image[index]){ min = image[index]; }
  }

  return min;
} // getMin()

/*
  The following method gets the maximum pixel value of the inputted image
*/
int getMax(unsigned char image[], int width, int height)
{
  int max = image[0];
  int index;
  for(index = 1; index < height*width; index++){
    if(max < image[index]){ max = image[index]; }
  }

  return max;
} // getMax()

/*
  The following method computes an histogram of the image, pixel by pixel
*/
int* buildHistogram(unsigned char image[], int width, int height, int min, int max)
{
  // Allocate memory for the histogram and initialise it to an array of 0s
  int* histogram = (int*) malloc(max * sizeof(int));
  int index;
  for(index = min; index <= max; index++){
    histogram[index] = 0;
  }

  // For every pixel value in the image
  // increment the correspondent column in the histogram
  int currentPixelValue;
  for(index = 0; index < height*width; index++){
    currentPixelValue = image[index];
    histogram[currentPixelValue]++;
  }

  return histogram;
} // buildHistogram()


/*
  The following method find the threshold from the histogram
  Looping from the minimum value to the maximum one
*/
int findThreshold(int* histogram, int min, int max)
{
  int currentThreshold;
  for (currentThreshold = min; currentThreshold <= max; currentThreshold++){

    // initialise variables
    int valueLess = 0;
    int valueGreater = 0;
    int counterLess = 0;
    int counterGreater = 0;
    double miuLess = 0;
    double miuGreater = 0;
    double average;

    // For all the values in the histogram below and equal to the threshold
    int histogramBrightness;
    for (histogramBrightness = min; histogramBrightness <= currentThreshold; histogramBrightness++){
      // Add the number of pixels that have that value of the histogram to the pixel counter
      counterLess += histogram[histogramBrightness];
      // Add the value of pixels that have that value of the histogram to the value counter
      valueLess += (histogram[histogramBrightness] * histogramBrightness);
    }

    // For all the values in the histogram above the threshold
    for (histogramBrightness = currentThreshold + 1; histogramBrightness <= max; histogramBrightness++){
      // Add the number of pixels that have that value of the histogram to the pixel counter
      counterGreater += histogram[histogramBrightness];
      // Add the value of pixels that have that value of the histogram to the value counter
      valueGreater += (histogram[histogramBrightness] * histogramBrightness);
    }

    // If there is at least one pixel below the threshold, compute the average below the threshold (miu l)
    if (counterLess != 0) { miuLess  = (double)valueLess  / (double)counterLess; }
    // If there is at least one pixel above the threshold, compute the average below the threshold (miu g)
    if (counterGreater != 0) { miuGreater  = (double)valueGreater  / (double)counterGreater; }

    // Average the two averages
    average = (miuLess + miuGreater) / 2;
    // Round them, since the threshold is an int
    average = round(average);

    // If the rounded average and the threshold are the same, then we have found the threshold
    if (average == currentThreshold) { return currentThreshold; }
  }
  return -1;
}// findThreshold()

/*
  The following method applies a given threshold to a given image
  The result is a binary image
*/
unsigned char* imageAfterThreshold(unsigned char* image, int width, int height, int threshold)
{
  // For any pixel in the image
  int currentPixel;
  for(currentPixel = 0; currentPixel < height*width; currentPixel++){
    // If the pixel is below or equal to the threshold, set it to 0
    if (image[currentPixel] <= threshold) { image[currentPixel] = 0;}
    // If the pixel is above the threshold, set it to 255
    else { image[currentPixel] = 255;}
  }
  return image;
}// imageAfterThreshold()

/*
  The following method compares two given ints
*/
int compare_ints (const void *a, const void *b)
{
  const int *ia = (const int *) a;
  const int *ib = (const int *) b;
  return (*ia > *ib) - (*ia < *ib);
}// compare_ints()

/*
  The following method applies median filtering to a given image.
  As suggested by Tim Morris, I am working with a 3x3 window.
  Median filtering is applied only on pixels that can be at the centre of the
  3x3 window (i.e. it can't be applied on the angles).
*/
unsigned char* medianFiltering(unsigned char* image, int width, int height)
{
  // Allocate memory for the final image and initialise it as the initial one
  unsigned char* smoothedImage = (unsigned char*) malloc((width * height) * sizeof(unsigned char));
  int currentPixel;
  for (currentPixel = 0; currentPixel < width * height; currentPixel++){
    smoothedImage[currentPixel] = image[currentPixel];
  }

  // I am working with a 3x3 window and I allocate its memory
  int windowSize = 3 * 3;
  int* window = (int*) malloc(windowSize * sizeof(int));
  // The current pixel is at the centre of the window, thus the edge is 1
  int edge = 1;

  int x;
  int y;
  // For every x and y that can be within a window
  for(y = edge; y < height - edge; y++){
      for(x = edge; x < width - edge; x++){
        // Compute the pixel
        currentPixel = (y * width) + x;
        // Get the current pixel
        window[0]= image[currentPixel];
        // Get the one on its left
        window[1]= image[currentPixel - edge];
        // Get the one on its right
        window[2]= image[currentPixel + edge];
        // Get the one on the left of the above one
        window[3]= image[currentPixel - width - edge];
        // Get the one above
        window[4]= image[currentPixel - width];
        // Get the one on the right of the above one
        window[5]= image[currentPixel - width + edge];
        // Get the one on the left of the below one
        window[6]= image[currentPixel + width - edge];
        // Get the one below
        window[7]= image[currentPixel + width];
        // Get the one on the right of the below one
        window[8]= image[currentPixel + width + edge];
        // Sort them
        qsort(window, windowSize, sizeof(int), compare_ints);
        // the element in the middle becomes the new value of the pixel
        smoothedImage[currentPixel] = window[4];
      }
  }
  // Deallocate memory of the window
  free(window);

  return smoothedImage;
} // medianFiltering()

/*
  The following method inverts the threshold values
  If the background is white, it turns black
  and viceversa
*/
void invertThreshold(unsigned char* image, int size){
    int index;
    for(index = 0; index < size; index++){
      if(image[index] == 0){ image[index] = 255;}
      else{ image[index] = 0; }
    }
}

/*
  The following method copies an array into a new one
*/
void copyArray(unsigned char* image_in, unsigned char* image_out, int size){
    int index;
    for(index = 0; index < size; index++){
      image_out[index] = image_in[index];
    }
}

/*
  The following method performs the equivalences of classes
  If the neighbour is not part of the background
  but the neighbour and the current pixel have different labels
  they are part of the same component and the two labels must belong to the same class
*/
void updateEquivalence(unsigned char* class, int* singleLabel, unsigned char neighbour_label, unsigned char current_label, unsigned char newLabel){
    int BACKGROUND = 0;
    int label;
  
    // If the neighbour is not part of the background
    // but the neighbour and the current pixel have different labels
    if ((neighbour_label != BACKGROUND) && (class[current_label] != class[neighbour_label])){
      // If the neighbour label is a newly entered label (== newLabel)
      // So it does not belong to a component
          if (singleLabel[class[neighbour_label]]) {
        // Its class becomes the one of the current one
            class[neighbour_label] = class[current_label];
        // The label of the neighbour is no longer single - it belongs to a component
              singleLabel[class[current_label]] = 0;
          }
      // If the current label is a newly entered label (== newLabel)
      // So it does not belong to a component
        else if (singleLabel[class[current_label]]){
        // Its class becomes the one of its neighbour
            class[current_label] = class[neighbour_label];
        // The label of the current one is no longer single - it belongs to a component
              singleLabel[class[current_label]] = 0;
          }
      // If they both belong to two different components
      // The two components are merged in just one
          else {
        // Starting from 1, since 0 corresponds to the background
            for (label = 1; label <= newLabel; label++)
          // For any label that belongs to the same class of the neighbour
          // Thus also the neighbour itself
                if (class[label] == class[neighbour_label])
            // Assign to it the class of the current label
                    class[label] = class[current_label];
          }
      }
} // updateEquivalence()

/*
  The following method makes the equivalence table to be contiguous
*/
void reorderEquivalences(unsigned char* class, unsigned char newLabel){
    int *temp = (int *) malloc ((newLabel+1) * sizeof(int));
      int tempLabel = 0;
  
    int label;
      for (label = 1; label <= newLabel; label++) {
      // If a class is equal to the label
      // It has not been used during the algorithm
          if (class[label] == label) {
              temp[label] = ++tempLabel;
          }
      }
  
    // class can now be filled with contiguous values
    // which are stored in temp
      for (label = 1; label <= newLabel; label++)
          class[label] = temp[class[label]];
  }
  
  /*
    The following method performs a connected component labelling algorithm on a given image
    The background is labelled 0.
    The components are consecutively labelled 1,2,3,etc..
  */
  unsigned char* CCA(unsigned char* image, int width, int height){
    int sizeOfImage = width*height;
    //initialise data structure for final output
    unsigned char* image_cca = (unsigned char*) malloc(sizeOfImage * sizeof(unsigned char));
    //initialise data structure
    int* singleLabel = (int*) malloc(sizeOfImage * sizeof(int));
    // background must be zero - threshold must be inverted
    invertThreshold(image, sizeOfImage);
    copyArray(image, image_cca, sizeOfImage);
    int index;
    for(index = 0; index < sizeOfImage; index++){
      singleLabel[index] = 0;
    }
  
    // equivalences are processed directly, so that equivalence classes are always updated
    // class[i] is the equivalence class associated with label i
    unsigned char* class = (unsigned char*) malloc(sizeOfImage * sizeof(unsigned char));
    // Initially esch label is assumed to belong to a different class
    int label;
    for(label = 0; label < (sizeOfImage); label++){
      class[label] = label;
    }
  
    unsigned char BACKGROUND = 0;
    unsigned char FOREGROUND = 255;
  
    unsigned char newLabel = 0;
  
    // p q r
    // s x
    unsigned char label_north_west = 0;    // p
    unsigned char label_north = 0;         // q
    unsigned char label_north_east = 0;    // r
    unsigned char label_west = 0;          // s
    unsigned char label_current_pixel = 0; // x
  
    int x_coord;
    int y_coord;
  
    //FIRST SCAN
  
    for(y_coord=0; y_coord < height; y_coord++){
      for(x_coord=0; x_coord < width; x_coord++){
  
        label_north_west = 0;    // p
        label_north = 0;         // q
        label_north_east = 0;    // r
        label_west = 0;          // s
        label_current_pixel = 0; //
  
        // If it is not background
        // Look for neighbours
        if (image[(y_coord * width) + x_coord] != BACKGROUND){
          // North-West
          if ((x_coord > 0) && (y_coord > 0))
            label_north_west = image_cca[((y_coord-1) * width) + x_coord-1];
          // North
          if (y_coord > 0)
            label_north = image_cca[((y_coord - 1) * width) + x_coord];
          // North-East
          if ((x_coord < width-1) && (y_coord > 0))
            label_north_east = image_cca[((y_coord-1) * width) + x_coord+1];
            // West
          if (x_coord > 0)
            label_west = image_cca[(y_coord * width) + x_coord-1];
  
          // All its neighbours are in the background - it's a new component
          if((label_north_west == BACKGROUND) && (label_north == BACKGROUND) && (label_north_east == BACKGROUND) && (label_west == BACKGROUND)){
            // Increment the label - get a new one
            if (newLabel == 255){ newLabel = 0;}
            newLabel++;
            // the new label is assigned to the current pixel
            label_current_pixel = newLabel;
            // The label does not belong to a component yet - it's single
            singleLabel[label_current_pixel] = 1;
            // The new label is temporarily set to be also the class of the current pixel
            // It will be merged with a component afterwards
            class[label_current_pixel] = label_current_pixel;
          }
          // Neighbours are in the foreground
              else {
            // new label is p
            label_current_pixel = label_north_west;
  
            // p is in the background - new label is q
            if (label_current_pixel == BACKGROUND)
              label_current_pixel = label_north;
  
                  // p and q are in the background - new label is r
                  if (label_current_pixel == BACKGROUND)
              label_current_pixel = label_north_east;
            // p/q and r are part of the same component - equivalence found
            else
                      updateEquivalence(class, singleLabel, label_north_east, label_current_pixel, newLabel);
  
                    // p,q and r are in the background - new label is s
                    if (label_current_pixel == BACKGROUND)
              label_current_pixel = label_west;
            else
              updateEquivalence(class, singleLabel, label_west, label_current_pixel, newLabel);
                  }
          image_cca[(y_coord*width)+x_coord] = label_current_pixel;
              }
          }
      } // FIRST SCAN
  
    // Reorder equivalence table to make it contiguous
    reorderEquivalences(class, newLabel);
  
      //SECOND SCAN
    //Substitue each label with the class it belongs to
      for (int y_coord = 0; y_coord < height; y_coord++) {
          for (int x_coord = 0; x_coord < width; x_coord++) {
              if (image_cca[(y_coord*width)+x_coord] > 0)
                  image_cca[(y_coord*width)+x_coord] = class[image_cca[(y_coord*width)+x_coord]];
          }
      }
  
    printf("\n The number of components in this image is: %d\n", getMax(image_cca, width, height));
    return image_cca;
} // CCA()

unsigned char* erosion(unsigned char* image, int width, int height)
{
  // Allocate memory for the final image and initialise it as the initial one
  unsigned char* erodedImage = (unsigned char*) malloc((width * height) * sizeof(unsigned char));
  int currentPixel;
  // I am working with a 3x3 window and I allocate its memory
  int windowSize = 3 * 3;
  int* window = (int*) malloc(windowSize * sizeof(int));
  // The current pixel is at the centre of the window, thus the edge is 1
  int edge = 1;

  int x;
  int y;
  // For every x and y that can be within a window
  for(y = edge; y < height - edge; y++){
    for(x = edge; x < width - edge; x++){
      // Compute the pixel
      currentPixel = (y * width) + x;
      // If the current pixel is not in the background
      if (image[currentPixel] == 0){
        // Get the one on its left
        window[0]= image[currentPixel - edge];
        // Get the one on its right
        window[1]= image[currentPixel + edge];
        // Get the one on the left of the above one
        window[2]= image[currentPixel - width - edge];
        // Get the one above
        window[3]= image[currentPixel - width];
        // Get the one on the right of the above one
        window[4]= image[currentPixel - width + edge];
        // Get the one on the left of the below one
        window[5]= image[currentPixel + width - edge];
        // Get the one below
        window[6]= image[currentPixel + width];
        // Get the one on the right of the below one
        window[7]= image[currentPixel + width + edge];
        // If at leats one neighbour is 0, set the pixel to 0
        int index;
        for(index = 0; index < 8; index++){
          if(window[index] == 255){
            erodedImage[currentPixel] = 255;
            break;
          }
        }
      }
      else {
        erodedImage[currentPixel] = 255;
      }
    }
  }
  // Deallocate memory of the window
  free(window);
  return erodedImage;
}

unsigned char* subtract(unsigned char* image, unsigned char* eroded, int width, int height)
{
  // Allocate memory for the final image and initialise it as the initial one
  unsigned char* finalImage = (unsigned char*) malloc((width * height) * sizeof(unsigned char));
  int currentPixel;
  for (currentPixel = 0; currentPixel < width * height; currentPixel++){
    finalImage[currentPixel] = eroded[currentPixel] - image[currentPixel];
  }
  return finalImage;
}

unsigned char* new_CCA(unsigned char* image, int width, int height)
{
  unsigned char* image_cca = (unsigned char*) malloc((width * height) * sizeof(unsigned char));
  int* neighbours = (int*) malloc((3*((width*height)+1)) * sizeof(int));

  int newLabel = 0;
  int y;
  int x;

  for (y = 0; y < height; y++){
    for (x = 0; x < width; x++){
      if ((image[(y*width)+x] != 0) && (image_cca[(y*width)+x] == 0)){
        newLabel++;
        label(image, image_cca, width, height, x, y, neighbours, newLabel);
      }
    }
  }
  return image_cca;
} //CCA()

#define x_coord (neighbours[stack_pointer-2])
#define y_coord (neighbours[stack_pointer-1])

#define assignLabel(x,y,label){         \
  neighbours[stack_pointer] = label;    \
  neighbours[stack_pointer+1] = x;      \
  neighbours[stack_pointer+2] = y;      \
  stack_pointer += 3;                   \
  goto recursive_routine;               \
}

#define propagateLabel {                      \
  stack_pointer -= 3;                         \
  int neighbour = neighbours[stack_pointer];  \
  switch (neighbour){                         \
    case 1 : goto east;                       \
    case 2 : goto north;                      \
    case 3 : goto south;                      \
    case 4 : goto done;                       \
    default: return;                          \
  }                                           \
}

void label(unsigned char* image, unsigned char* image_cca, int width, int height, int x, int y, int* neighbours, int newLabel)
{
  // LABEL
  neighbours[0] = 0;
  // X COORDINATE
  neighbours[1] = x;
  // Y COORDINATE
  neighbours[2] = y;

  int stack_pointer = 3;
  int currentPixel;

  recursive_routine:
    currentPixel = (y_coord*width)+x_coord;
    // pixel is in the background - propagate
    if (image [currentPixel] == 0)
      propagateLabel;
    // pixel has lareayd been labelled - propagate
    if (image_cca[currentPixel] != 0)
      propagateLabel;

    image_cca[currentPixel] = newLabel;
    // WEST
    if (x_coord > 0)
      assignLabel(x_coord-1, y_coord, 1);

    // EAST
    east:
    if (x_coord < width-1)
      assignLabel(x_coord+1, y_coord, 2);

    // NORTH
    north:
    if (y_coord > 0)
      assignLabel(x_coord, y_coord-1, 3);

    // SOUTH
    south:
    if (y_coord < height-1)
      assignLabel(x_coord, y_coord+1, 4);

    done:
  propagateLabel;
} // label()

int main(int argc, char *argv[]){
    if (argc < 8) {
        printf("Usage: ./image_processing sickle.jpg sickle_thresholded.jpg sickle_smoothed.jpg sickle_cca.jpg sickle_eroded.jpg sickle_perimeter.jpg sickle_perimeter_cca.jpg\n");
        return 1;
    }

    unsigned char *image;
    int width, height, channels;
    read_JPEG_file(argv[1], &width, &height, &channels, &image);

    // get min value in image
    int minValue = getMin(image, width, height);
    // get max value in image
    int maxValue = getMax(image, width, height);
    // compute histogram of image
    int* histogram = buildHistogram(image, width, height, minValue, maxValue);
    // find best threshold
    int threshold = findThreshold(histogram, minValue, maxValue);
    printf("\n The computed threshold for this image is: %d\n\n", threshold);

    // Apply thresholding function
    image = imageAfterThreshold(image, width, height, threshold);

    // Write thresholded image to output file
    write_JPEG_file(argv[2], width, height, channels, image, 95);

    // Apply median filtering to thresholded image
    image_median = medianFiltering(image, width, height);

    // Write smoothed image to output file
    write_JPEG_file(argv[3], width, height, channels, image_median, 95);

    // Apply CCA algorithm
    image = CCA(image_median, width, height);

    // Write labelled image to output file
    write_JPEG_file(argv[4], width, height, channelabel_west, image, 95);

    // Apply erosion to thresholded image
    unsigned char* erodedimage = erosion(image_median, width, height);

    // Write eroded image to output file
    write_JPEG_file(argv[5], width, height, channels, erodedimage, 95);

    // Image with perimeters
    unsigned char* perimeterimage = subtract(image_median, erodedimage, width, height);

    // Write perimeter image to output file
    write_JPEG_file(argv[6], width, height, channels, perimeterimage, 95);

    unsigned char* labelledperimeterimage = CCA(perimeterimage, width, height);

    // Write perimeter image to output file
    write_JPEG_file(argv[7], width, height, channels, labelledperimeterimage, 95);

    free(image);
    free(histogram);

    return 0;
}

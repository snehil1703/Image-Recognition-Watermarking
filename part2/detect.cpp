//
// detect.cpp : Detect cars in satellite images.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h> 
#include <iostream>

using namespace std;

#define PI 3.14159265

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &cars)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &cars, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  SDoublePlane output(input.rows(), input.cols());
  
  // Convolution code here
  
  if (filter.cols() != filter.rows() || filter.cols() % 2 != 1) {
    cout << "Bad dimensions for filter in function convolve_general";
  }

  long i, j, u, v, k;
  long use_i, use_j;
  double sum;

  k = (filter.cols() - 1) / 2;

  for (i = 0; i < input.rows(); i++) {
    for (j = 0; j < input.cols(); j++) {
      sum = 0;
      for (u = -k; u <= k; u++) {
        for (v = -k; v <= k; v++) {
          // handle image edges
          // assume pixel values past the edge of the image are equal to their nearest pixel...

          
          if ((i - u) < 0) {
            use_i = 0; 
          } else if ((i - u) >= input.rows()) {
            use_i = input.rows() - 1;
          } else {
            use_i = i - u;
          }

          if ((j - v) < 0) {
            use_j = 0;
          } else if ((j - v) >= input.cols()) {
            use_j = input.cols() - 1;
          } else {
            use_j = j - v;
          }
          
          
          sum += input[use_i][use_j] * filter[k+u][k+v];
        }
      }
      output[i][j] = sum;
    }
  }

  return output;
}


// Convolve an image with a separable convolution kernel
// //
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane output2(input.rows(), input.cols());

  long i, j, u, v, k;
  long use_i, use_j;
  double sum;

  k = (row_filter.cols() - 1) / 2;

  for (i = 0; i < input.rows(); i++) {
    for (j = 0; j < input.cols(); j++) {
      sum = 0;
      
      for (v = -k; v <= k; v++) {

        if ((j - v) < 0) {
          use_j = 0;
        } else if ((j - v) >= input.cols()) {
          use_j = input.cols() - 1;
        } else {
          use_j = j - v;
        }

        sum += input[i][use_j] * row_filter[0][k+v];
      }
      
      output[i][j] = sum;
    }
  }


  for (i = 0; i < input.rows(); i++) {
    for (j = 0; j < input.cols(); j++) {
      sum = 0;
      for (u = -k; u <= k; u++) {
        if ((i - u) < 0) {
          use_i = 0;
        } else if ((i - u) >= input.rows()) {
          use_i = input.rows() - 1;
        } else {
          use_i = i - u;
        }
        sum += output[use_i][j] * col_filter[k+u][0];
 
      
      }
      output2[i][j] = sum;
    }
  }

  return output2;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane output2(input.rows(), input.cols());
  // Implement a sobel gradient estimation filter with 1-d filters
 
  //set up row and column vectors for vertical edge detection with separable convolution.
  SDoublePlane col_filter(3,1);
  col_filter[0][0] = 1;
  col_filter[1][0] = 2;
  col_filter[2][0] = 1;

  SDoublePlane row_filter(1,3);
  row_filter[0][0] = -1;
  row_filter[0][1] = 0;
  row_filter[0][2] = 1;

  //run the convolution
  output = convolve_separable(input, row_filter, col_filter); 

  //set up vectors for horizontal edge detection
  col_filter[0][0] = -1;
  col_filter[1][0] = 0;
  col_filter[2][0] = 1;

  row_filter[0][0] = 1;
  row_filter[0][1] = 2;
  row_filter[0][2] = 1;

  //run the convolution again
  output2 = convolve_separable(input, row_filter, col_filter); 

  
  for (int i=0; i < output.rows(); i++) {
    for (int j=0; j < output.cols(); j++) {
      output[i][j] = pow( pow(output[i][j], 2) + pow(output2[i][j], 2), .5);
    }
  }
  
  return output;
}

SDoublePlane clip_image(const SDoublePlane &input, int row, int col, int height, int width) 
{
  SDoublePlane output(height, width);

  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      output[i][j] = input[row+i][col+j];
    }
  }
  
  return output;
}

SDoublePlane threshold(const SDoublePlane &input, double top_thresh, double second_thresh) 
{
  SDoublePlane output(input.rows(), input.cols());
 
  for (int i = 0; i < input.rows(); i++) {
    for (int j = 0; j < input.cols(); j++) {
       if (input[i][j] >= top_thresh) {
        output[i][j] = 255;
      } else if (input[i][j] >= second_thresh) {
        output[i][j] = 0;
      } else {
        output[i][j] = 0;
      }
    }
  }
  return output;
}

double get_gradient(double x1, double x2, double y1, double y2) 
{
  if (x2 - x1 == 0) {
    return 0;
  }

  return atan((y2 - y1) / (x2 - x1));
}

SDoublePlane nonmax_suppression(const SDoublePlane &input) 
{
  SDoublePlane output(input.rows(), input.cols());
  double grad;
  double val1, val2;
 
  for (int i = 1; i < input.rows() - 1; i++) {
    for (int j = 1; j < input.cols() - 1; j++) {
      grad = get_gradient(input[i+1][j+1], input[i][j], input[i][j+1], input[i+1][j]);

      if (grad <= M_PI / 8 || grad > 7 * M_PI / 8) {
        val1 = input[i+1][j+1];
        val2 = input[i-1][j-1]; 
      } else if (grad <= 3 * M_PI / 8) {
        val1 = input[i][j+1];
        val2 = input[i][j-1];
      } else if (grad <= 5 * M_PI / 8) {
        val1 = input[i-1][j+1];
        val2 = input[i+1][j-1];
      } else /* (grad <= 7 * M_PI / 8) */ {
        val1 = input[i-1][j];
        val2 = input[i+1][j];
      }

      if (input[i][j] >= val1 && input[i][j] >= val2) {
        output[i][j] = input[i][j];
      } else {
        output[i][j] = 0;
      }
    }
  }

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  
  return output;
}


int distanceToOrigin(int i, int j,int angle_deg) {
  double angle = M_PI * angle_deg / 180;

  double x1 = j;
  double y1 = i;
  double x2 = j + cos(angle);
  double y2 = i + sin(angle);
  double x3 = 0;
  double y3 = 0;

  double v1x = x2 - x1;
  double v1y = y2 - y1;
  double v2x = x3 - x1;
  double v2y = y3 - y1;

  double dot_product = v1x * v2x + v1y * v2y;

  double px = x1 + dot_product * v1x;
  double py = y1 + dot_product * v1y;

  double x_intercept;
  int intercept_mult;

  if (angle == 0) {
    intercept_mult = 1;
  } else {
    x_intercept = x1 - v1x * (y1 / v1y);
    //cout << "x intercept is " << x_intercept;
    if (x_intercept >= 0) {
      intercept_mult = 1;
    } else {
      intercept_mult = -1;
    }
  }
  
  return (int) intercept_mult * pow(px * px + py * py,.5); 
}

SDoublePlane hough(const SDoublePlane &input)
{

  int max_angle = 180;
  SDoublePlane output(2*(input.rows() + input.cols()), max_angle);
  int zero_dist_row = input.rows() + input.cols();

  for (int i = 1; i < input.rows() - 1; i++) {
    for (int j = 1; j < input.cols() - 1; j++) {
      if (input[i][j] > 0) {
           int angle;
           for (angle=0;angle<max_angle;angle++){
             int dist=distanceToOrigin(-i,j,angle);
             output[zero_dist_row + dist][angle] += 1;
           } 	
      }
    }
  }

  return output;
}


void get_hough_col_info(const SDoublePlane &input) {
  //for hough transform, sum of every column should be the same...
  double col_sum = 0;

  int max_threshold = (input.rows() + input.cols()) / 40;
  
  for (int i = 0; i < input.rows(); i++) {
    col_sum += input[i][0];
  }
  
  //cout << "Column sum is " << col_sum << "\n";

  double ave = col_sum / input.rows();
 
  double all_sd_sum = 0;
  
  for (int j = 0; j < input.cols(); j++) {
    double sd_sum = 0;
    int num_linies = 0;
    int num_maxes = 0;
    double sum_diffs = 0;
    for (int i = 0; i < input.rows(); i++) {
      sd_sum += pow(input[i][j] - ave, 2);
      if (input[i][j] >= 27) {
        num_linies++;
      }
    }

    for (int i = 1; i < input.rows() - 1; i++) {
      if (input[i][j] > input[i-1][j] && input[i][j] >= input[i+1][j] && input[i][j] > max_threshold) {
        num_maxes++;
      }
      sum_diffs += abs(input[i][j] - input[i+1][j]);
    }

    double sd = pow(sd_sum / input.rows(), .5);
    //cout << "SD of col " << j << " is " << sd << "\n";
    //cout << "Linies for col " << j << " is " << num_linies << "\n";
    //cout << "Maxes for col " << j << " is " << num_maxes << "\n";
    //cout << "Sum of diffs " << j << " is " << sum_diffs << "\n";
    all_sd_sum += sd;
  }
  //cout << "Sum of all sds : " << all_sd_sum << "\n";
}

double get_highest_diff_ratio(const SDoublePlane &input) {
  double highest_diff = 0;
  double lowest_diff = 100000;

  for (int j = 0; j < input.cols(); j++) {
    double sum_diffs = 0;
    for (int i = 0; i < input.rows() - 1; i++) {
      sum_diffs += abs(input[i][j] - input[i+1][j]);
    }
    if (sum_diffs > highest_diff) {
      highest_diff = sum_diffs;
    }
    if (sum_diffs < lowest_diff) {
      lowest_diff = sum_diffs;
    } 
  }
  
  return (highest_diff / lowest_diff);
}

double get_orth_ratio(const SDoublePlane &input, int angle) {
  double first_diff = 0;
  double second_diff = 0;

  int second_angle = (angle < 90) ? angle + 90 : angle - 90;

  for (int i = 0; i < input.rows() - 1; i++) {
    first_diff += abs(input[i][angle] - input[i+1][angle]);
    second_diff += abs(input[i][second_angle] - input[i+1][second_angle]);
  }
  if (second_diff == 0) {
    second_diff = .0001;
  }
  return (first_diff / second_diff);
}

//get the rows and values of local maxima in a column (of a hough trans.)
void get_maxes_for_col(const SDoublePlane &input, int the_col) {
  int j = the_col;
  int diff = (input.cols() + input.rows()) / 2;

  for (int i = 1; i < (input.rows() - 1); i++) {
    if (input[i][j] > input[i-1][j] && input[i][j] >= input[i+1][j] ) {
      //cout << "Col " << the_col << " has max at row " << i << " of " << input[i][j] << "\n";
    }
  }
}

int get_hough_max_col(const SDoublePlane &input) {
  double best_val = 0;
  int best_col = 0;

  for (int j = 0; j < input.cols(); j++) {
    for (int i = 0; i < input.rows(); i++) {
      if (input[i][j] > best_val) {
         best_val = input[i][j];
         best_col = j;
      } 
    }
  }
  //cout << "Column with highest val ... " << best_col << "!!!!\n";
  return best_col;
}

long count_edge_pixels(const SDoublePlane &input) {
  long sum_pixels = 0;
  for (int j = 0; j < input.cols(); j++) {
    for (int i = 0; i < input.rows(); i++) {
      if (input[i][j] > 0) {
         sum_pixels++;
      }
    }
  }
  return sum_pixels;

}

//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  
  // test step 2 by applying mean filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;

  SDoublePlane vert_blur_filter(3,3);
  vert_blur_filter[0][1] = 1/3.0;
  vert_blur_filter[1][1] = 1/3.0;
  vert_blur_filter[2][1] = 1/3.0;
  
 
  //SDoublePlane output_image = convolve_general(input_image, mean_filter); 
  SDoublePlane output_image = sobel_gradient_filter(input_image, true); 
  
  output_image = nonmax_suppression(output_image);

  output_image = threshold(output_image, 150, 55);

  SDoublePlane hough_plane = hough(output_image);


  //multiple box blur
  for (int i; i < 10; i++) {
    hough_plane =  convolve_general(hough_plane, vert_blur_filter);
  }

  // randomly generate some detected cars -- you'll want to replace this
  //  with your car detection code obviously!
  vector<DetectedBox> cars;

  int win_size = 30;
  int max_i = (input_image.rows() - win_size) / win_size;
  int max_j = (input_image.cols() - win_size) / win_size;

  for (int i = 0; i <= max_i; i++) {
    for (int j = 0; j <= max_j; j++) {
      int row = i * win_size;
      int col = j * win_size;

      SDoublePlane input_clipped = clip_image(output_image, row, col, win_size, win_size);
      SDoublePlane hough_clipped = hough(input_clipped);
     
      int best_col = get_hough_max_col(hough_clipped);

      double ratio = get_highest_diff_ratio(hough_clipped);

      double orth = get_orth_ratio(hough_clipped, best_col);

      double orth_to_low = ratio / orth;
      
      long edge_pixels = count_edge_pixels(input_clipped);
      double log_regression = -.0203 * ratio + .013 * orth + -.00256 * orth_to_low + .00279 * edge_pixels;
  
      //cout << ratio << ", " << orth << ", " << orth_to_low << ", " << edge_pixels << ", " << log_regression << "\n";

      if (log_regression > .35) {
        DetectedBox s;
        s.row = row;
        s.col = col;
        s.width = win_size;
        s.height = win_size;
        s.confidence = 1;
        cars.push_back(s);
      }
    }
  }

  write_detection_txt("detected.txt", cars);
  write_detection_image("detected.png", cars, input_image);
  write_detection_image("hough.png", cars, hough_plane);
  write_detection_image("edges.png", cars, output_image);
}


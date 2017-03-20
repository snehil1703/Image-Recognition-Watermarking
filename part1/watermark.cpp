//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// Snehil Vishwakarma
//
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>
#include <SImage.h>
#include <SImageIO.h>
#include <fft.h>

using namespace std;

#define loop(i, n) for( ; i < n ; ++i )
#define loop0(i, n) for( i = 0 ; i < n ; ++i )
#define same3(t) t, t, t
#define sub(a, b) (a - b)
#define add(a, b) (a + b)
#define sqr(a) (a * a)

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
void add_padding(SDoublePlane &input)
{   
    int i = input.rows(); 
    if (i < input.cols())
        i = input.cols();
    int rcmax = 1;
    while (rcmax < i)
        rcmax <<= 1;
    if ( rcmax > input.rows() || rcmax > input.cols() )
    {
      double val = 0, valt;
      int r = input.rows(), c = input.cols(), ct=0;
      double rc = (double)r*(double)c;
      for (int x = 0; x<r; x++)
      {
          for(int y =0; y<c;y++)
          {
              val += input[x][y];
              valt = fmod( val , rc);
              if (val != valt)
                  ct++;
              val = valt;
          }
      }
      while (val > 0)
          val/=10.0;
      valt = (double)ct + val;
      input.paddingval(rcmax,val);
    }
}

// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
    fft_real = input;
    fft_imag = SDoublePlane(input.rows(), input.cols());

    FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
    output_real = input_real;
    SDoublePlane output_imag = input_imag;

    FFT_2D(0, output_real, output_imag);
}

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag)
{
    int i, j, max = fft_real.rows();
    SDoublePlane specto( max, max);
  
    int mn_x, mn_y, mx_x, mx_y;
    loop0 (i,max)
        loop0 (j, max)
        {  
            specto[i][j] = log(sqrt(add(sqr(fft_real[i][j]),sqr(fft_imag[i][j])))); 
            if(i == 0 && j == 0)
            {
                mn_x = i;
                mn_y = j;
                mx_x = i;
                mx_y = j;
            }
            if(specto[mn_x][mn_y] > specto[i][j])
            {
                mn_x = i;
                mn_y = j;
            }
            if(specto[mx_x][mx_y] < specto[i][j])
            {
                mx_x = i;
                mx_y = j;
            }
        }
    double mx = specto[mx_x][mx_y], mn = specto[mn_x][mn_y];
    loop0 (i,max)
        loop0 (j,max)
            specto[i][j] = (specto[i][j]-mn) * 255 / (mx-mn) ;
    return specto;
}

// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input)
{
    int i, j, k, max = input.rows();
    bool chk = false;
    SDoublePlane fft_r, fft_i;
    fft(input, fft_r, fft_i);
    SDoublePlane clean( max, max);
  
    int mn_x, mn_y, mx_x, mx_y;
    loop0 (i, max)
    {
        loop0 (j, max)
        {  
            clean[i][j] = log(sqrt(add(sqr(fft_r[i][j]),sqr(fft_i[i][j]))));
            if(i == 0 && j == 0)
            {
                mn_x = i;
                mn_y = j;
                mx_x = i;
                mx_y = j;
            }
            if(clean[mn_x][mn_y] > clean[i][j])
            {
                mn_x = i;
                mn_y = j;
            }
            if(clean[mx_x][mx_y] < clean[i][j])
            {
                mx_x = i;
                mx_y = j;
            }
        }
    }
  
    double mx = clean[mx_x][mx_y], mn = clean[mn_x][mn_y];
    loop0 (i,max)
        loop0 (j,max)
            clean[i][j] = (clean[i][j]-mn) * 255 / (mx-mn) ;
    
    loop0 (i,max)
    {  
        loop0 (j, max)
        {
      
            if(clean[i][j] > 128 && !(i > (max/2)-80 && i < (max/2)+80) && !(j > (max/2)-80 && j < (max/2)+80))
            {
                //cout << i<< " " << j << " " << clean[i][j] << endl;
                k = 1;
                int m, n, x, y, z;
                chk = false;
                while(1)
                { 
                    m = i+k;
                    n = j+k;
                    x = i-k;
                    z = j-k;
                    while(x < 0)
                        x++;
                    while(z < 0)
                        z++;
                    while(m > max)
                        m--;
                    while(n > max)
                        n--;
                    if(x == 0 && z == 0 && m == max && n == max)
                        break;
                    loop ( x, m)
                    {
                        y = z;
                        loop ( y, n)
                        {
                            if(clean[x][y] <= 128)
                            {
                                chk = true;
                                break;
                            }
                        }
                        if(chk)
                            break;
                    }
                    if(chk)
                        break;
                    k++;
                }
                //cout << i << " " << j << " new: " << x <<" "<< y << endl;
                clean[i][j] = clean[x][y];
                fft_r[i][j] = fft_r[x][y];
                fft_i[i][j] = fft_i[x][y];
            }     
        }        
    }

    clean = SDoublePlane( max, max);
    ifft( fft_r, fft_i, clean);
    SDoublePlane clean1 = clean;
    //double sum = 0;
    loop0 (i,max)
        loop0 (j, max)
        {  
	          //sum += clean[i][j];
	          //sum = fmod(sum,256.0);

            if(clean[i][j]-128 >= 0)
		            clean[i][j] = clean[i][j]*10;
	          else
		            clean[i][j] = clean[i][j]/10;
        }
    //cout << endl << sum << endl;
    //sum /= 128;
    //sum = (double)(259 * (4+255)) / (double)(255 * (259-4));
    //cout << endl << sum << endl;
    /*loop0 (i,max)
    {
        loop0 (j,max)
        {
	          clean[i][j] = sum*(clean[i][j]-128)+128;
	          if(clean[i][j] < 0)
		            clean[i][j] = 0;
	          if(clean[i][j] > 255)
		            clean[i][j] = 255;
        }
    }*/
    SDoublePlane specto = fft_magnitude( fft_r, fft_i);
    string output_specto = "pictures/noise1/spectogram_cleaned_noise1.png";
    SImageIO::write_png_file((output_specto).c_str(), same3(specto));
    
    string output_opt = "pictures/noise1/opt_cleaned_noise1.png";
    SImageIO::write_png_file((output_opt).c_str(), same3(clean));
    cout << "Optimized Out: " << output_opt << endl;
    
    return clean1;
}

//Random Values Generator
vector<bool> rng(int N, int l)
{
    srand(N);
    vector<bool> ans;
    while(l > 0)
    {
        ans.push_back( (bool)(rand()%2));
        l--;
    }
    return ans;
}

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N)
{
    int d = input.rows();
    
    SDoublePlane fft_r, fft_i;
    fft(input, fft_r, fft_i);
    int l = d/2;
    //if (l*l != d)
    //    l = (int)sqrt(d/2);
    
    double min,max;
    min = fft_r[0][0];
    max = fft_r[0][0];
    for (int i=0;i<d;i++)
    {
        for(int j=0;j<d;j++)
        {
            if(max < fft_r[i][j])
                max = fft_r[i][j];
            if(min > fft_r[i][j])
                min = fft_r[i][j];
        }
    }
    double alpha = 2.5;
    
    vector<bool> bv = rng(N, l);
    
    //double r = (double)l/4;
    double c = (double)(l-1)/2;
    
    double temp = ((double)l/2) - 0.5;
    
    int mx = (int)(c-0.5);
    int my = (int)(c-temp);
    int px = (int)(c+0.5);
    int py = (int)(c+temp);
    for( int step = 0; step < (l/2); step++, mx--, my++, px++, py--)
    {
	      //cout << mx << " " << my << endl; 
        //cout << fft_r[mx][my] << " " << fft_r[my][px] << endl;
        //cout << std::abs(fft_r[mx][my]) << " " << std::abs(fft_r[mx][my]) << endl;
        //cout << (double)bv[step] << " " << (double)bv[step + (l/2)] << endl;
        fft_r[mx][my] = fft_r[mx][my] + (alpha * std::abs(fft_r[mx][my]) * (double)bv[step]);
        fft_r[my][px] = fft_r[my][px] + (alpha * std::abs(fft_r[my][px]) * (double)bv[step + (l/2)]);
        fft_r[px][py] = fft_r[mx][my];
        fft_r[py][mx] = fft_r[my][px];
        //cout << fft_r[mx][my] << " " << fft_r[my][px] << endl;
    }
    
    //SDoublePlane specto = fft_magnitude( fft_r, fft_i);
    //string output = "spectogram_watermark.png";
    //SImageIO::write_png_file((output).c_str(), same3(specto));
    
    SDoublePlane answer( d, d);
    ifft( fft_r, fft_i, answer);
    return answer;
}

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N)
{
    int d = input.rows();
    SDoublePlane fft_r, fft_i;
    fft(input, fft_r, fft_i);
    int l = d/2;
    //if (l*l != d)
    //  l = (int)sqrt(d/2);
    
    vector<bool> bv = rng(N, l);
    
    double c = (double)(l-1)/2;
  
    double temp = (((double)l)/2) - 0.5;
  
    int mx = (int)(c-0.5);
    int my = (int)(c-temp);
    int px = (int)(c+0.5);
    int py = (int)(c+temp);
    double xbar = 0, ybar = 0;
    for( int step = 0; step < (l/2); step++, mx--, my++, px++, py--)
    {
        xbar += ( fft_r[mx][my] + fft_r[my][px] );
        ybar += ( bv[step] + bv[step+(l/2)] );
    }
    xbar = xbar / l;
    ybar = ybar / l;
    mx = (int)(c-0.5);
    my = (int)(c-temp);
    px = (int)(c+0.5);
    py = (int)(c+temp);
    double xy = 0, x2 = 0, y2 = 0;
    for( int step = 0; step < (l/2); step++, mx--, my++, px++, py--)
    {
        //cout << mx << " " << my << endl; 
        //cout << fft_r[mx][my] << " " << fft_r[my][px] << endl;

        xy += ( (fft_r[mx][my] - xbar) * (bv[step] - ybar) + (fft_r[my][px] - xbar) * (bv[step+(l/2)] - ybar) );
        x2 += ( sqr(fft_r[mx][my] - xbar) + sqr(fft_r[my][px] - xbar) );
        y2 += ( sqr(bv[step] - ybar) + sqr(bv[step+(l/2)] - ybar) );
    }
    double r = xy / sqrt( x2 * y2 );
    
    //cout << endl << r << " ";
    
    SDoublePlane answer;
    if ( r >= 0.6)
    {
        //cout << "Watermark present. YES! \n";
        answer = SDoublePlane( d, d);
        ifft( fft_r, fft_i, answer);
    }else
    {    
        //cout << "Watermark absent. NO! \n";
        answer = SDoublePlane( d-1, d-1);
    }
    return answer;
}

int main(int argc, char **argv)
{
    try 
    {
        if(argc < 4)
        {
            cout << "Insufficent number of arguments; correct usage:" << endl;
	          cout << "    p2 problemID inputfile outputfile" << endl;
	          return -1;
        }
    
        string part = argv[1];
        string inputFile = argv[2];
        string outputFile = argv[3];
        outputFile = inputFile.substr( 0, inputFile.find_last_of("/")) + "/" + outputFile;
        
        SDoublePlane input_image = SImageIO::read_png_file(inputFile.c_str());
        add_padding(input_image);
    
        if(part == "1.1")
        {
            SDoublePlane fft_r, fft_i;
            fft(input_image, fft_r, fft_i);
            SDoublePlane specto = fft_magnitude( fft_r, fft_i);
            SImageIO::write_png_file(outputFile.c_str(), same3(specto));
            cout << "In: " << inputFile <<"  Out: " << outputFile << endl << endl;
        }
        else if(part == "1.2")
        {
            cout << "In: " << inputFile <<"  Out: " << outputFile << endl << endl;
            SDoublePlane clean = remove_interference(input_image);
            SImageIO::write_png_file(outputFile.c_str(), same3(clean));
        }
        else if(part == "1.3")
        {
            if(argc < 6)
	          {
	              cout << "Need 6 parameters for watermark part:" << endl;
	              cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	              return -1;
	          }
	          string op(argv[4]);
	          int N = atoi(argv[5]);
	          if(op == "add")
	          {
	              SDoublePlane marked = mark_image(input_image, N);
	              //cout << marked.rows() << " " << marked.cols() << " " << input_image.rows() << " " << input_image.cols() << endl;
	              SImageIO::write_png_file(outputFile.c_str(), same3(marked));
	              cout << "In: " << inputFile <<"  Out: " << outputFile << endl << endl;
	          }
	          else if(op == "check")
	          {
	              SDoublePlane checked = check_image(input_image, N);
	              string temp = inputFile.substr( inputFile.find_last_of("/"));
	              temp = temp.substr( temp.find_last_of("_")+1);
	              temp = temp.substr( 0, temp.find_last_of("."));
	              int trueN = atoi(temp.c_str());
	              if (checked.rows() == input_image.rows())
	              {
	                  if (trueN != N)
                        cout << " !!False Positive!! True Value: " << trueN << " Accepted Value: " << N << endl;
	                  else
	                      cout << "Watermark present. YES! " << N << endl;
	                  
	              }
	              else
	              {   
	                  if (trueN == N)
	                      cout << "!Watermark missed! " << N << endl;
	              }
	              //SImageIO::write_png_file(outputFile.c_str(), same3(checked));
	              //cout << "In: " << inputFile <<"  Out: " << outputFile << endl;
	          }
	          else
	              throw string("Bad operation!");
        }
        else
            throw string("Bad part!");
    } 
    catch(const string &err) 
    {
        cerr << "Error: " << err << endl;
    }
}

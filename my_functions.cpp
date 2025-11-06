#include "my_functions.hpp"
#include <math.h>
#include <numeric>  // Necessário para std::accumulate
#include <vector>   // Necessário para std::vector

// Define constants for the multi-level histogram
#define COARSE_BINS 16 // 256 / 16
#define FINE_BINS 16   // 4-bit coarse, 4-bit fine


#include <omp.h>
// if OpenMP is used, the command to compile must have -fopenmp

void my_sobel_parallelized (Mat im_in, Mat im_out)
{
    int rows = im_in.rows;
    int cols = im_in.cols;

    #pragma omp parallel for shared(im_in, im_out, rows, cols)

    for (int i = 1; i< rows-1; i++){
		for (int j = 1; j< cols-1; j++){
			int n = im_in.at<uchar>(i-1, j);
			int s = im_in.at<uchar>(i+1, j);
			int e = im_in.at<uchar>(i, j+1);
			int w = im_in.at<uchar>(i,j-1);
			int ne = im_in.at<uchar>(i-1,j+1);
			int nw = im_in.at<uchar>(i-1,j-1);
			int se = im_in.at<uchar>(i+1,j+1);
			int sw = im_in.at<uchar>(i+1,j-1);
			int c = im_in.at<uchar>(i,j);

        	int gx = 2*e + ne + se - 2*w - sw - nw;
			int gy = 2*n + ne + nw - 2*s - sw - se;
			int g = abs(gx) + abs(gy);
			if (g>255)
				g=255;
			im_out.at<uchar>(i, j) = (uchar) (g);
        }
    }
}

void my_sobel(Mat im_in, Mat im_out)
{

 int rows = im_in.rows; // height
 int cols = im_in.cols; // width
 
 
    for (int i = 1; i< cols-1; i++)
 	for (int j = 1; j< rows-1; j++){
 		int n = im_in.at<uchar>(j-1,i);
 		int s = im_in.at<uchar>(j+1,i);
 		int e = im_in.at<uchar>(j,i+1);
 		int w = im_in.at<uchar>(j,i-1);
 		int ne = im_in.at<uchar>(j-1,i+1);
 		int nw = im_in.at<uchar>(j-1,i-1);
 		int se = im_in.at<uchar>(j+1,i+1);
 		int sw = im_in.at<uchar>(j+1,i-1);
 		int c = im_in.at<uchar>(j,i);
 	
 	 	int gx = 2*e + ne + se - 2*w - sw - nw + 0*n + 0*c + 0*s;
 	 	int gy = 2*n + ne + nw - 2*s - sw - se + 0*w + 0*c + 0*e;
 	 	
 	 	float g = sqrt ((float)(gx*gx) + (float)(gy*gy));
 	 	
 	 	if (g>255)
 	 		g=255;
 	 	
 	 	im_out.at<uchar>(j,i) = (uchar) (g);
 	 	
 	 	}

};


void my_sobel_optimized(const Mat& im_in, Mat& im_out)
{
	const int rows = im_in.rows; // height
	const int cols = im_in.cols; // width
	const uchar* __restrict lastLine = im_in.ptr<uchar>(0);
    const uchar* __restrict curLine = im_in.ptr<uchar>(1); 
    uchar* __restrict outLine;
	const uchar* __restrict nxtLine;

    for (int i = 1; i < rows - 1; i++)
    {
        nxtLine = im_in.ptr<uchar>(i + 1);
        outLine = im_out.ptr<uchar>(i);

		int nw = lastLine[0];
        int w = curLine[0]; 
        int sw = nxtLine[0]; 

        int n = lastLine[1];
        int c = curLine[1];
        int s = nxtLine[1];

        __builtin_prefetch(nxtLine + 64, 0, 1);

        for (int j = 1; j < cols - 1; j++)
        {
            int ne = lastLine[j + 1]; // Coluna j+1 (norte)
            int e = curLine[j + 1];  // Coluna j+1 (centro)
            int se = nxtLine[j + 1];  // Coluna j+1 (sul)z'

            int gx = 2*e + ne + se - 2*w - sw - nw;
            int gy = 2*n + ne + nw - 2*s - sw - se;
            int g = abs(gx) + abs(gy);
            if (g > 255)
                g = 255;
                
            outLine[j] = (uchar)(g);

			// move the sliding window
			nw = n;
			w = c;
			sw = s;
			n = ne;
			c = e;
			s = se;
        }
        lastLine = curLine;
        curLine = nxtLine;
    }
};

void my_median (Mat im_in, Mat im_out, int n)
{

 int k = n*n;
 int rows = im_in.rows; // height
 int cols = im_in.cols; // width
 int r, c, rr, cc, p;
 vector<uchar> v (k);
 
    for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         p = 0;
         for(rr=(r-(n/2));rr<(r-(n/2)+n);rr++){
            for(cc=(c-(n/2));cc<(c-(n/2)+n);cc++){
               if((rr>=0)&&(rr<rows)&&(cc>=0)&&(cc<cols)){
                  v[p] = im_in.at<uchar>(rr,cc);
                  p++;
	       }
            }
         }
        
         sort (v.begin(), v.end());
         im_out.at<uchar>(r,c) = v[k/2+1];
      }

}

};

void my_median_optimized(const Mat& im_in, Mat& im_out, int n)
{
    const int r = n / 2;
    const int rows = im_in.rows, cols = im_in.cols;
    const int cols_pad = cols + 2*r;

    CV_Assert(im_in.type() == CV_8UC1);
    im_out.create(rows, cols, CV_8UC1);

    // 1) Pad once
    cv::Mat im_pad;
    copyMakeBorder(im_in, im_pad, r, r, r, r, cv::BORDER_REFLECT_101);

    // 2) Column histograms (flattened): Hcol[j*256 + v]
    std::vector<uint16_t> Hcol(cols_pad * 256, 0);
    auto Hcol_at = [&](int j)->uint16_t* { return &Hcol[j * 256]; };

    // 3) Initialize first n rows into each column histogram
    for (int y = 0; y < n; ++y) {
        const uchar* __restrict p = im_pad.ptr<uchar>(y);
        for (int j = 0; j < cols_pad; ++j) {
            Hcol_at(j)[p[j]]++;
        }
    }

    // 4) Kernel histogram
    uint16_t Hker[256];

    const int median_thr = (n*n + 1) / 2;
    for (int y = 0; y < rows; ++y) {
        // slide column histograms one row down (except y==0)
        if (y > 0) {
            const uchar* __restrict top = im_pad.ptr<uchar>(y - 1);
            const uchar* __restrict bot = im_pad.ptr<uchar>(y + n - 1);
            for (int j = 0; j < cols_pad; ++j) {
                uint16_t* __restrict hc = Hcol_at(j);
                hc[ top[j] ]--;
                hc[ bot[j] ]++;
            }
        }

        // build kernel histogram for j=0: sum of first n columns
        std::fill_n(Hker, 256, 0);
        for (int j = 0; j < n; ++j) {
            const uint16_t* __restrict hc = Hcol_at(j);
            #pragma GCC unroll 256
            for (int v = 0; v < 256; ++v) Hker[v] += hc[v];
        }

        uchar* __restrict out = im_out.ptr<uchar>(y);

        // horizontal sliding
        for (int x = 0; x < cols; ++x) {
            // find median
            int acc = 0, v = 0;
            #pragma GCC unroll 256
            for (; v < 256; ++v) { acc += Hker[v]; if (acc >= median_thr) break; }
            out[x] = (uchar)v;

            if (x == cols - 1) break;

            // Hker = Hker - Hcol[x] + Hcol[x+n]
            const uint16_t* __restrict hL = Hcol_at(x);
            const uint16_t* __restrict hR = Hcol_at(x + n);
            #pragma GCC unroll 256
            for (int t = 0; t < 256; ++t) Hker[t] = Hker[t] - hL[t] + hR[t];
        }
    }
}
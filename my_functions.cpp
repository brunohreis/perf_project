#include "my_functions.hpp"
#include <math.h>


void my_sobel (Mat im_in, Mat im_out)
{
	int rows = im_in.rows; // height
	int cols = im_in.cols; // width
	uchar* lastLine = im_in.ptr<uchar>(0);
    uchar* curLine = im_in.ptr<uchar>(1); 
    uchar* outLine;
	uchar* nxtLine;

    for (int i = 1; i < rows - 1; i++)
    {
        nxtLine = im_in.ptr<uchar>(i + 1);
        outLine = im_out.ptr<uchar>(i);
        for (int j = 1; j < cols - 1; j++)
        {
            int nw = lastLine[j - 1];
			int n  = lastLine[j];
            int ne = lastLine[j + 1];

            int w  = curLine[j - 1];
            int e  = curLine[j + 1];

            int sw = nxtLine[j - 1];
            int s  = nxtLine[j];
            int se = nxtLine[j + 1];

            int gx = 2*e + ne + se - 2*w - sw - nw;
            int gy = 2*n + ne + nw - 2*s - sw - se;
            int g = abs(gx) + abs(gy);
            if (g > 255)
                g = 255;
                
            outLine[j] = (uchar)(g);
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
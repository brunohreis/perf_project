#include "my_functions.hpp"
#include <math.h>
#include <numeric>  // Necessário para std::accumulate
#include <vector>   // Necessário para std::vector

// Define constants for the multi-level histogram
#define COARSE_BINS 16 // 256 / 16
#define FINE_BINS 16   // 4-bit coarse, 4-bit fine


/* #include <omp.h>
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
} */


void my_sobel(const Mat& im_in, Mat& im_out)
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

/* void my_median_better(Mat im_in, Mat im_out, int n){

	int r = n/2;
	int rows = im_in.rows; // height
	int cols = im_in.cols; // width
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

	for(int i=0; i<rows; i++){
		vector<vector<int>> histograms;
		// inicialização de um vetor de histogramas de colunas, de dimensão [2r+1][2r+1]
		// preenchimento dos histogramas (primeiras 2r+1 colunas, começando pela linha i)
		// cálculo e atribuição em im_out da mediana do elemento im_out[r+1][r+1] utilizando histogramas e bucket_sort O(1)
		for(int j=1; j<cols; j++){
			// computa-se o histograma da coluna a direita do kernel, adição do elemento inferior a última coluna do kernel, e subtração do elemento superior (no histograma da última coluna) 
			// subtração do primeiro histograma do vetor de histogramas
			// computação da mediana em O(1) e atribuição em im_out[r+1][r+1]
		}
	}
}; */

// --- 1. NEW DATA STRUCTURES ---

/**
 * @brief Holds a 2-level histogram (Section III-C)
 * 'coarse' has 16 bins (for upper 4 bits)
 * 'fine' has 256 bins (for all 8 bits)
 */
struct MultiLevelHist {
    vector<int> coarse; // 16 bins [cite: 185]
    vector<int> fine;   // 256 bins [cite: 184]

    MultiLevelHist() : coarse(COARSE_BINS, 0), fine(256, 0) {}
};

// --- 2. NEW HELPER FUNCTIONS (for MultiLevelHist) ---

/**
 * @brief Adds a pixel value to a multi-level histogram
 */
inline void hist_add_multi(MultiLevelHist& hist, uchar value) {
    hist.coarse[value >> 4]++; // Add to coarse bin (upper 4 bits) [cite: 183-185]
    hist.fine[value]++;        // Add to fine bin (all 8 bits)
}

/**
 * @brief Subtracts a pixel value from a multi-level histogram
 */
inline void hist_sub_multi(MultiLevelHist& hist, uchar value) {
    hist.coarse[value >> 4]--; // Sub from coarse bin
    hist.fine[value]--;        // Sub from fine bin
}

/**
 * @brief Adds only the COARSE level of a column histogram to the kernel's coarse histogram
 */
inline void hist_add_coarse(vector<int>& H_kernel_coarse, const MultiLevelHist& h_col) {
    for (int k = 0; k < COARSE_BINS; k++) {
        H_kernel_coarse[k] += h_col.coarse[k];
    }
}

/**
 * @brief Subtracts only the COARSE level
 */
inline void hist_sub_coarse(vector<int>& H_kernel_coarse, const MultiLevelHist& h_col) {
    for (int k = 0; k < COARSE_BINS; k++) {
        H_kernel_coarse[k] -= h_col.coarse[k];
    }
}

inline uchar find_median_on_demand(
    const vector<int>& H_kernel_coarse, // 16 bins, always up-to-date
    vector<int>& H_kernel_fine_segment, // 16 bins, REBUILT every time
    const vector<MultiLevelHist>& col_hists,
    int median_threshold,
    int j_local, // Current *local* column index (relative to stripe)
    int n        // Kernel width
) {
    // 1. Find median in COARSE histogram [cite: 187-189]
    int sum_coarse = 0;
    int coarse_idx = 0;
    for (; coarse_idx < COARSE_BINS; coarse_idx++) {
        sum_coarse += H_kernel_coarse[coarse_idx];
        if (sum_coarse >= median_threshold) {
            break;
        }
    }

    // 2. Calculate local threshold for the fine segment
    int fine_threshold = median_threshold - (sum_coarse - H_kernel_coarse[coarse_idx]);
    int fine_segment_start_bin = coarse_idx * FINE_BINS;

    // 3. ALWAYS Rebuild the 16-bin fine segment from scratch
    // This fixes the visual bug by never using stale data.
    H_kernel_fine_segment.assign(FINE_BINS, 0); // Reset buffer to all zeros
    
    int start_col = j_local; // Left-most column in the window *within the stripe buffer*
    int end_col = j_local + n;
    
    for (int k_col = start_col; k_col < end_col; k_col++) {
        for (int k_bin = 0; k_bin < FINE_BINS; k_bin++) {
            H_kernel_fine_segment[k_bin] += col_hists[k_col].fine[fine_segment_start_bin + k_bin];
        }
    }

    // 4. Find median in the now-correct FINE segment
    int sum_fine = 0;
    int fine_idx = 0;
    for (; fine_idx < FINE_BINS; fine_idx++) {
        sum_fine += H_kernel_fine_segment[fine_idx];
        if (sum_fine >= fine_threshold) {
            break;
        }
    }

    // 5. Return the final median value
    return (uchar)(fine_segment_start_bin + fine_idx);
}


// --- 4. FUNÇÃO PRINCIPAL DA MEDIANA (OTIMIZADA PARA CACHE) ---

// Define a largura do bloco para ser "amigável" ao cache L2
// Um histograma tem ~1KB. 128 hist. ~= 128KB, que cabe facilmente no cache.
#define CACHE_FRIENDLY_STRIPE_WIDTH 128

void my_median_better(Mat im_in, Mat im_out, int n)
{
    int r = n / 2; 
    int k_size = n * n; 
    int median_threshold = (k_size / 2) + 1;
    int rows = im_in.rows;
    int cols = im_in.cols;

    // 1. Crie a imagem com preenchimento (padding) UMA VEZ
    Mat im_pad;
    copyMakeBorder(im_in, im_pad, r, r, r, r, BORDER_REFLECT_101);

    // 2. Otimização (Seção III-B): Loop sobre "Blocos Verticais"
    // Isso garante que o vetor 'col_hists' caiba no cache 
    for (int j_base = 0; j_base < cols; j_base += CACHE_FRIENDLY_STRIPE_WIDTH) 
    {
        int j_start = j_base;
        int j_end = std::min(j_base + CACHE_FRIENDLY_STRIPE_WIDTH, cols);
        int stripe_width = j_end - j_start;

        // 3. Inicialize um VETOR DE HISTOGRAMA PEQUENO (que cabe no cache)
        // Precisamos de 'stripe_width' colunas + 'n' colunas extras para a janela
        vector<MultiLevelHist> col_hists(stripe_width + n);

        // Loop principal (linha por linha) *dentro* do loop de bloco
        for (int i = 0; i < rows; i++) {
            
            if (i == 0) {
                // Primeira linha: Preenche os histogramas do bloco do zero
                for (int j_s = 0; j_s < (stripe_width + n); j_s++) {
                    int j_pad = j_start + j_s; // Coluna global na imagem com padding
                    if (j_pad >= im_pad.cols) break;
                    
                    col_hists[j_s] = MultiLevelHist(); // Reseta o histograma
                    for (int i_k = 0; i_k < n; i_k++) {
                        hist_add_multi(col_hists[j_s], im_pad.at<uchar>(i_k, j_pad));
                    }
                }
            } else {
                // Outras linhas: "Desliza" os histogramas do bloco para baixo
                for (int j_s = 0; j_s < (stripe_width + n); j_s++) {
                    int j_pad = j_start + j_s;
                    if (j_pad >= im_pad.cols) break;
                    
                    uchar val_remove = im_pad.at<uchar>(i - 1, j_pad);
                    uchar val_add    = im_pad.at<uchar>(i + n - 1, j_pad);
                    hist_sub_multi(col_hists[j_s], val_remove);
                    hist_add_multi(col_hists[j_s], val_add);
                }
            }

            // 5. Inicializa o histograma 'coarse' do kernel para a primeira coluna do bloco
            vector<int> kernel_hist_coarse(COARSE_BINS, 0);
            for (int j_k = 0; j_k < n; j_k++) {
                hist_add_coarse(kernel_hist_coarse, col_hists[j_k]);
            }

            // Buffer temporário para o segmento fino (agora dentro do loop de linha)
            vector<int> kernel_hist_fine_segment(FINE_BINS, 0);

            // 6. Loop de Coluna (processa apenas colunas *deste* bloco)
            for (int j = j_start; j < j_end; j++) {
                
                // j_local é o índice no 'col_hists' (de 0 a stripe_width)
                int j_local = j - j_start; 
                
                // 7. Calcula a mediana (usando a função CORRIGIDA)
                im_out.at<uchar>(i, j) = find_median_on_demand(
                    kernel_hist_coarse, kernel_hist_fine_segment, col_hists, 
                    median_threshold, j_local, n
                );

                if (j == j_end - 1) break; // Não desliza na última coluna do bloco
                
                // 8. "Desliza" o histograma 'coarse' do kernel
                hist_sub_coarse(kernel_hist_coarse, col_hists[j_local]);
                hist_add_coarse(kernel_hist_coarse, col_hists[j_local + n]);
            }
        }
    }
}
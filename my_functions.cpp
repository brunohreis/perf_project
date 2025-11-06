#include "my_functions.hpp"
#include <math.h>
#include <numeric> 
#include <vector>  

void my_sobel(Mat im_in, Mat im_out)
{
	const int rows = im_in.rows; // height
	const int cols = im_in.cols; // width
	const uchar* __restrict lastLine = im_in.ptr<uchar>(0);
    const uchar* __restrict curLine = im_in.ptr<uchar>(1); 
    /*uchar* outLine;
	uchar* nxtLine;*/

    for (int i = 1; i < rows - 1; ++i)
    {
        const uchar* __restrict nxtLine = im_in.ptr<uchar>(i + 1);
        uchar* __restrict outLine = im_out.ptr<uchar>(i);

		int nw = lastLine[0];
        int w = curLine[0]; 
        int sw = nxtLine[0]; 

        int n = lastLine[1];
        int c = curLine[1];
        int s = nxtLine[1];

        // prefetch for the future line
        __builtin_prefetch(nxtLine + 64, 0, 1);

        for (int j = 1; j < cols - 1; ++j)
        {
            int ne = lastLine[j + 1]; // Coluna j+1 (norte)
            int e = curLine[j + 1];  // Coluna j+1 (centro)
            int se = nxtLine[j + 1];  // Coluna j+1 (sul)z'

            int gx = 2*e + ne + se - 2*w - sw - nw;
            int gy = 2*n + ne + nw - 2*s - sw - se;
            int g = abs(gx) + abs(gy);
            /*
            if (g > 255)
                g = 255;
            */
            outLine[j] = (uchar)(g > 255 ? 255 : g);

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

/**
 * @brief Adiciona um valor de pixel a um histograma (vetor de 256 posições).
 */
inline void hist_add(vector<int>& hist, uchar value) {
    hist[value]++;
}

/**
 * @brief Subtrai um valor de pixel de um histograma.
 */
inline void hist_sub(vector<int>& hist, uchar value) {
    hist[value]--;
}

/**
 * @brief Adiciona um histograma (h_add) a outro (H_kernel).
 * Isto é o O(256) mencionado no artigo[cite: 150].
 */
inline void hist_add_hist(vector<int>& H_kernel, const vector<int>& h_add) {
    for (int k = 0; k < 256; k++) {
        H_kernel[k] += h_add[k];
    }
}

/**
 * @brief Subtrai um histograma (h_sub) de outro (H_kernel).
 * Isto é o O(256) mencionado no artigo[cite: 151].
 */
inline void hist_sub_hist(vector<int>& H_kernel, const vector<int>& h_sub) {
    for (int k = 0; k < 256; k++) {
        H_kernel[k] -= h_sub[k];
    }
}

/**
 * @brief Encontra o valor mediano em um histograma.
 * Esta é a operação O(1) (na verdade O(128) em média)[cite: 85, 152].
 * @param hist O histograma do kernel.
 * @param threshold O valor para encontrar ( (n*n) / 2 + 1 ).
 * @return O valor do pixel mediano (uchar).
 */
inline uchar find_median_from_hist(const vector<int>& hist, int threshold) {
    int sum = 0;
    for (int k = 0; k < 256; k++) {
        sum += hist[k];
        if (sum >= threshold) {
            return (uchar)k;
        }
    }
    return 255; // Caso não encontre (não deve acontecer)
}


// ======================================================================
// FUNÇÃO PRINCIPAL MY_MEDIAN (O(1) com Histogramas de Coluna)
// ======================================================================

void my_median_better(Mat im_in, Mat im_out, int n)
{
    int r = n / 2; // Raio do kernel
    int k_size = n * n; // Tamanho total do kernel (ex: 9, 25)
    int median_threshold = (k_size / 2) + 1; // Posição da mediana

    int rows = im_in.rows;
    int cols = im_in.cols;

    // 1. Criar uma imagem 'im_pad' com bordas replicadas.
    // Isso elimina *todos* os `if` de verificação de borda.
    Mat im_pad;
    copyMakeBorder(im_in, im_pad, r, r, r, r, BORDER_REFLECT_101);

    // 2. Inicialização do vetor de histogramas de colunas 
    // Um histograma (de 256 posições) para CADA coluna da imagem PADDED.
    vector<vector<int>> col_hists(im_pad.cols, vector<int>(256, 0));

    // 3. Histograma do Kernel (H)
    vector<int> kernel_hist(256, 0);

    // Loop principal (linha por linha)
    for (int i = 0; i < rows; i++) {
        
        // ============================================================
        // A. INICIALIZAÇÃO DA LINHA 'i'
        // ============================================================
        
        // 4. Preenchimento/Atualização dos histogramas de coluna
        if (i == 0) {
            // -- Primeira linha (i=0): Preenchimento inicial --
            // Preenche todos os histogramas de coluna com as primeiras 'n' linhas
            for (int j_pad = 0; j_pad < im_pad.cols; j_pad++) {
                for (int i_k = 0; i_k < n; i_k++) {
                    hist_add(col_hists[j_pad], im_pad.at<uchar>(i_k, j_pad));
                }
            }
        } else {
            // -- Outras linhas (i>0): "Deslizar" todos os histogramas de coluna 1 linha para baixo --
            // Esta é a Etapa 1 (Figura 2a) do artigo 
            for (int j_pad = 0; j_pad < im_pad.cols; j_pad++) {
                // Remove o pixel de cima (que saiu da janela)
                uchar val_remove = im_pad.at<uchar>(i - 1, j_pad);
                // Adiciona o pixel de baixo (que entrou na janela)
                uchar val_add    = im_pad.at<uchar>(i + n - 1, j_pad);
                
                hist_sub(col_hists[j_pad], val_remove);
                hist_add(col_hists[j_pad], val_add);
            }
        }

        // 5. Cálculo do histograma do kernel (H) para a *primeira coluna* (j=0)
        // O kernel H é a soma dos 'n' primeiros histogramas de coluna 
        kernel_hist.assign(256, 0); // Limpa o H
        for (int j_k = 0; j_k < n; j_k++) {
            hist_add_hist(kernel_hist, col_hists[j_k]);
        }

        // ============================================================
        // B. LOOP DA COLUNA 'j' (Deslizando o kernel H)
        // ============================================================
        for (int j = 0; j < cols; j++) {
            
            // 6. Computação da mediana em O(1) e atribuição [cite: 138, 145]
            // (O kernel_hist atual representa a janela centrada em (i, j))
            im_out.at<uchar>(i, j) = find_median_from_hist(kernel_hist, median_threshold);

            // 7. Preparar o kernel_hist para a *próxima* iteração (j+1)
            // (Não faz isso na última coluna)
            if (j == cols - 1) {
                break; 
            }
            
            // "Desliza" o kernel H uma coluna para a direita 
            // Esta é a Etapa 2 (Figura 2b) do artigo
            
            // Subtrai o histograma da coluna que saiu (esquerda)
            hist_sub_hist(kernel_hist, col_hists[j]);
            
            // Adiciona o histograma da coluna que entrou (direita)
            hist_add_hist(kernel_hist, col_hists[j + n]);
        }
    }
}

void my_median_better_v2(const Mat& im_in, Mat& im_out, int n)
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

    for (int y = 0; y < rows; ++y) {
        // slide column histograms one row down (except y==0)
        if (y > 0) {
            const uchar* __restrict top = im_pad.ptr<uchar>(y - 1);
            const uchar* __restrict bot = im_pad.ptr<uchar>(y + n - 1);
            for (int j = 0; j < cols_pad; ++j) {
                Hcol_at(j)[ top[j] ]--;
                Hcol_at(j)[ bot[j] ]++;
            }
        }

        // build kernel histogram for j=0: sum of first n columns
        std::fill(std::begin(Hker), std::end(Hker), 0);
        for (int j = 0; j < n; ++j) {
            const uint16_t* __restrict hc = Hcol_at(j);
            #pragma GCC unroll 256
            for (int v = 0; v < 256; ++v) Hker[v] += hc[v];
        }

        // slide horizontally
        uchar* __restrict out = im_out.ptr<uchar>(y);
        const int median_thr = (n*n + 1) / 2;

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
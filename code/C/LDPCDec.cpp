#include "LDPCDec.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef FAST_MATH
#include "fastexp.h"
#include "fastlog.h"
#endif

#ifdef _WIN32
#define log2(x) 			(log(x)/log((double)2))
#endif

/// \brief	defines the maximum absolute value of L-values used
///
#define MAXLLR	100

/// \brief	returns minimum of a and b
///
#define MIN(a,b)			(((a) < (b)) ? (a) : (b))

/// \brief	returns maximum of a and b
///
#define MAX(a,b)			(((a) > (b)) ? (a) : (b))

/// \brief	clips a at magnitude b
///
#define CLIP(a,b)			MAX(MIN(a,b), -b)

/// \brief	return sign of a
///
#define SIGN(a)				(((a)<0) ? (-1) : (1))

#define FABS(a)				(((a)<0) ? (-a) : (a))

/// \brief	Standard constructor.
///

CLDPCDec::CLDPCDec(const char *filename) {
    m_N = 0;
    m_M = 0;
    m_E = 0;
    m_g = 0;
    m_vardegree = NULL;
    m_chkdegree = NULL;
    m_interleaver = NULL;
    m_v2c = NULL;
    m_c2v = NULL;
    m_hard = NULL;
    m_Lsyn = NULL;

    m_clearmsg = true;
    m_method = DECODER_SPA;

    m_decloaded = false;

    m_corrveclen = 0;
    m_corrvec = NULL;
    
    readDecoder(filename);
    
}

/// \brief	Standard destructor.
///

CLDPCDec::~CLDPCDec(void) {
    freeMem();
}

/// \brief  Allocate memory.
///

void CLDPCDec::allocateMem() {
    freeMem();
    m_vardegree = (unsigned int*) malloc(m_N * sizeof (unsigned int));
    m_chkdegree = (unsigned int*) malloc(m_M * sizeof (unsigned int));
    m_interleaver = (unsigned int*) malloc(m_E * sizeof (unsigned int));
    m_v2c = (double*) malloc(m_E * sizeof (double));
    m_c2v = (double*) malloc(m_E * sizeof (double));
    // m_hard is used in the method encodeRA where it needs to be one element
    // larger than the number of edges. This is because the last edge of the
    // accumulator is ommited in the code definition.
    m_hard = (int*) malloc((m_E + 1) * sizeof (int));
    m_Lsyn = (double*) malloc(m_M * sizeof (double));
    m_decloaded = true;
}

/// \brief  Free memory.
///

void CLDPCDec::freeMem() {
    free(m_vardegree);
    free(m_chkdegree);
    free(m_interleaver);
    free(m_v2c);
    free(m_c2v);
    free(m_hard);
    free(m_corrvec);
    free(m_Lsyn);

    m_vardegree = NULL;
    m_chkdegree = NULL;
    m_interleaver = NULL;
    m_v2c = NULL;
    m_c2v = NULL;
    m_hard = NULL;
    m_decloaded = false;
    m_corrvec = NULL;
    m_corrveclen = 0;
    m_Lsyn = NULL;
}

/// \brief	Returns the length of the codeword.
///

unsigned int CLDPCDec::getN() {
    return m_N;
}

/// \brief	Returns the number of parity-checks.
///

unsigned int CLDPCDec::getM() {
    return m_M;
}

/// \brief	Returns the number of information bits.
///

unsigned int CLDPCDec::getK() {
    return m_N - m_M;
}

/// \brief	Returns the number of edges.
///

unsigned int CLDPCDec::getE() {
    return m_E;
}

/// \brief  Decoding of an LDPC code using defined method
/// \param	Lch			L-values received from the channel
/// \param	Lapp		a-posteriori L-values after decoding
/// \param	maxit		maximum number of iterations to perform
///
/// On termination the number of iterations performed is returned.

unsigned int CLDPCDec::decode(double *Lch, double *Lapp, unsigned int maxit) {
    switch (m_method) {
        case DECODER_SPA:
            return decodeSPA(Lch, Lapp, maxit);
        case DECODER_MSA:
            return decodeMSA(Lch, Lapp, maxit);
        default:
            return decodeSPA(Lch, Lapp, maxit);
    }
}

/// \brief  Decoding of an LDPC code using the sum-product algorithm.
/// \param	Lch			L-values received from the channel
/// \param	Lapp		a-posteriori L-values after decoding
/// \param	maxit		maximum number of iterations to perform
///
/// This is the implementation of the sum-product algorithm for decoding
/// of an LDPC code. The decoder performs up to maxit iterations or
/// terminates if all check nodes are satisfied.
/// On termination the number of iterations performed is returned.

unsigned int CLDPCDec::decodeSPA(double *Lch, double *Lapp, unsigned int maxit) {
    unsigned int m, n, e, d;
    unsigned int j, i, it;
    int parity;
    double x;
    double d1, d2, i1;
    bool valid;

    if (m_clearmsg) {
        // initialize messages (in LLR domain)
        for (e = 0; e < m_E; e++)
            m_c2v[e] = 0.0;
    }

    // perform iterations
    valid = false;
    for (it = 0; it < maxit; it++) {
        // varnode processing
        d = 0;
        for (n = 0; n < m_N; n++) {
            Lapp[n] = Lch[n];
            for (j = 0; j < m_vardegree[n]; j++)
                Lapp[n] += m_c2v[d + j];

            for (j = 0; j < m_vardegree[n]; j++) {
#ifdef FAST_MATH
                m_v2c[d + j] = fasterexp(-CLIP(Lapp[n] - m_c2v[d + j], MAXLLR));
#else
                m_v2c[d + j] = exp(-CLIP(Lapp[n] - m_c2v[d + j], MAXLLR));
#endif
                m_hard[d + j] = (Lapp[n] <= 0);
            }

            d += m_vardegree[n];
        }

        // chknode processing
        valid = true;
        d = 0;
        for (m = 0; m < m_M; m++) {
            parity = 0;
            for (j = 0; j < m_chkdegree[m]; j++)
                parity ^= m_hard[m_interleaver[d + j]];

            if (parity)
                valid = false;

            for (j = 0; j < m_chkdegree[m]; j++) {
                x = 0.0;

                for (i = 0; i < m_chkdegree[m]; i++)
                    if (j != i) {
                        i1 = m_v2c[m_interleaver[d + i]];
                        d1 = x + i1;
                        d2 = 1 + x*i1;
                        x = d1 / d2;
                    }
                //x = (x+m_v2c[m_interleaver[d+i]])/(1+x*m_v2c[m_interleaver[d+i]]);
#ifdef FAST_MATH
                m_c2v[m_interleaver[d + j]] = -fasterlog(x);
#else
                m_c2v[m_interleaver[d + j]] = -log(x);
#endif			
            }
#ifdef FAST_MATH
            m_Lsyn[m] = -fasterlog((x + m_v2c[m_interleaver[d + m_chkdegree[m] - 1]]) / (1 + x * m_v2c[m_interleaver[d + m_chkdegree[m] - 1]]));
#else
            m_Lsyn[m] = -log((x + m_v2c[m_interleaver[d + m_chkdegree[m] - 1]]) / (1 + x * m_v2c[m_interleaver[d + m_chkdegree[m] - 1]]));
#endif 
            d += m_chkdegree[m];
        }

        if (valid)
            break;

    }

    return it;
}

/// \brief  Decoding of an LDPC code using the min-sum algorithm.
/// \param	Lch			L-values received from the channel
/// \param	Lapp		a-posteriori L-values after decoding
/// \param	maxit		maximum number of iterations to perform
///
/// This is the implementation of the min-sum algorithm for decoding
/// of an LDPC code. The decoder performs up to maxit iterations or
/// terminates if all check nodes are satisfied.
/// On termination the number of iterations performed is returned.

unsigned int CLDPCDec::decodeMSA(double *Lch, double *Lapp, unsigned int maxit) {
    unsigned int m, n, e, d;
    unsigned int j, i, it;
    int parity;
    double x, s;
    //int s;
    bool valid;

    if (m_clearmsg) {
        // initialize messages (in LLR domain)
        for (e = 0; e < m_E; e++)
            m_c2v[e] = 0.0;
    }

    // perform iterations
    for (it = 0; it < maxit; it++) {
        // varnode processing
        d = 0;
        for (n = 0; n < m_N; n++) {
            Lapp[n] = Lch[n];
            for (j = 0; j < m_vardegree[n]; j++)
                Lapp[n] += m_c2v[d + j];

            for (j = 0; j < m_vardegree[n]; j++) {
                m_v2c[d + j] = CLIP(Lapp[n] - m_c2v[d + j], MAXLLR);
                m_hard[d + j] = (Lapp[n] <= 0);
            }

            d += m_vardegree[n];
        }

        // chknode processing
        valid = true;
        d = 0;
        for (m = 0; m < m_M; m++) {
            parity = 0;
            for (j = 0; j < m_chkdegree[m]; j++)
                parity ^= m_hard[m_interleaver[d + j]];

            if (parity)
                valid = false;

            for (j = 0; j < m_chkdegree[m]; j++) {
                x = MAXLLR;
                s = 1;

                for (i = 0; i < m_chkdegree[m]; i++)
                    if (j != i) {
                        //x		= MIN(x, fabs(m_v2c[m_interleaver[d+i]]));
                        x = MIN(x, FABS(m_v2c[m_interleaver[d + i]]));
                        s *= SIGN(m_v2c[m_interleaver[d + i]]);
                    }

                m_c2v[m_interleaver[d + j]] = s*x;
                //m_c2v[m_interleaver[d+j]] = copysign(x,s);
            }

            //m_Lsyn[m] = MIN(x, fabs(m_v2c[m_interleaver[d+j]])) * s * SIGN(m_v2c[m_interleaver[d+j]]);
            m_Lsyn[m] = MIN(x, FABS(m_v2c[m_interleaver[d + j]])) * s * SIGN(m_v2c[m_interleaver[d + j]]);

            d += m_chkdegree[m];
        }

        // post-processing
        if (it < m_corrveclen) {
            for (e = 0; e < m_E; e++)
                m_c2v[e] *= m_corrvec[it];
        }

        if (valid)
            break;
    }

    return it;
}

/// \brief  Reads parity-check matrix from file.
/// \param	filename	name of the file
///
/// Reads a parity-check matrix from a file which is in our format.

int CLDPCDec::readDecoder(const char *filename) {
    FILE *file;
    unsigned int n, m, e;

    file = fopen(filename, "r");
    if (file == NULL)
        return 0;

    fscanf(file, "%u %u %u", &m_N, &m_M, &m_E);

    allocateMem();

    for (n = 0; n < m_N; n++)
        fscanf(file, "%u", &(m_vardegree[n]));
    for (m = 0; m < m_M; m++)
        fscanf(file, "%u", &(m_chkdegree[m]));
    for (e = 0; e < m_E; e++)
        fscanf(file, "%u", &(m_interleaver[e]));

    fclose(file);

    return 1;
}

/// \brief  Sets parameter clearmsg.
/// \param	clearmsg	new value
///
/// Sets parameter clearmsg.

void CLDPCDec::set_param_clearmsg(bool clearmsg) {
    m_clearmsg = clearmsg;
}

/// \brief  Sets parameter method.
/// \param	method	new value
///
/// Sets parameter method.

void CLDPCDec::set_param_method(unsigned int method) {
    m_method = method;
}

/// \brief  Sets correction vector for MSA post-processing
/// \param	len      number of correction vectors
/// \param  corrvec  pointer to vector
///
/// Sets parameter method.

void CLDPCDec::set_param_corrvec(unsigned int len, double *corrvec) {
    unsigned int i;

    free(m_corrvec);
    m_corrvec = (double*) malloc(len * sizeof (double));
    m_corrveclen = len;

    for (i = 0; i < m_corrveclen; i++)
        m_corrvec[i] = corrvec[i];
}

unsigned int CLDPCDec::encodeLDGM(unsigned int *info, unsigned int *parity) {
    unsigned int m, k, d;
    unsigned int j;

    // varnode processing
    d = 0;
    for (k = 0; k < getK(); k++) {
        for (j = 0; j < m_vardegree[k]; j++)
            m_hard[d + j] = info[k];

        d += m_vardegree[k];
    }

    // interleaving and chknode processing
    d = 0;
    for (m = 0; m < m_M; m++) {
        parity[m] = 0;
        for (j = 0; j < (m_chkdegree[m] - 1); j++)
            parity[m] ^= m_hard[m_interleaver[d + j]];

        d += m_chkdegree[m];
    }

    return 1;
}

unsigned int CLDPCDec::encodeRA(unsigned int *info, unsigned int *parity) {
    unsigned int m, k, dv, dc;
    unsigned int j;

    // varnode processing
    dv = 0;
    for (k = 0; k < getK(); k++) {
        for (j = 0; j < m_vardegree[k]; j++)
            m_hard[dv + j] = info[k];

        dv += m_vardegree[k];
    }

    // interleaving, chknode processing and accumulating
    dc = 0;
    for (m = 0; m < m_M; m++) {
        parity[m] = 0;
        for (j = 0; j < (m_chkdegree[m] - 1); j++)
            parity[m] ^= m_hard[m_interleaver[dc + j]];

        m_hard[dv + 1] = parity[m];

        dv += 2;
        dc += m_chkdegree[m];
    }

    return 1;
}

/// \brief  Estimates the mutual information from a vector of L-values
///
/// Estimates the mutual information from a vector of L-values of length N

double CLDPCDec::L2mutual(double *L, long N) {
    long n;
    double I = 1;

    for (n = 0; n < N; n++)
        I -= log2(1 + exp(-L[n])) / N;

    return I;
}

/// \brief  Mutual information of the syndrome
///
/// Estimates the mutual information of the syndrome

double CLDPCDec::syndromeInformation() {
    return L2mutual(m_Lsyn, m_M);
}

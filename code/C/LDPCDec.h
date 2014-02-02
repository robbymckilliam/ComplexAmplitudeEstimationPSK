#ifndef LDPCDEC_H
#define LDPCDEC_H

#define DECODER_SPA		0
#define DECODER_MSA		1

#include <string>

/// \brief		Main class for LDPC decoding
/// \author		Gottfried Lechner (gottfried.lechner@unisa.edu.au)
/// \version	3.0
/// \date     Created      : March 2008
/// \date     Last modified: April 2012
///
/// This class holds the LDPC decoder.

class CLDPCDec {
public:
    CLDPCDec(const char* filename);
    ~CLDPCDec(void);

    int readDecoder(const char* filename);

    unsigned int encodeLDGM(unsigned int *info, unsigned int *parity);
    unsigned int encodeRA(unsigned int *info, unsigned int *parity);
    unsigned int decode(double *Lch, double *Lapp, unsigned int maxit);

    void set_param_clearmsg(bool clearmsg);
    void set_param_method(unsigned int method);
    void set_param_corrvec(unsigned int len, double *corrvec);

    double syndromeInformation();

    /// Determines whether the decoder is loaded
    bool dec_loaded() {
        return m_decloaded;
    }

    unsigned int getN() const;
    unsigned int getM() const;
    unsigned int getK() const;
    unsigned int getE() const;

    ///Compute the log likelihood ration for binary phase shift keying with input x and
    ///given amplitude and noise variance
    static double BPSK2LLR(double x, double var);
    
    ///Return expected BPSK value given a log likelihood ratio
    static double LLR2BPSK(double llr);

private:
    unsigned int m_N;
    unsigned int m_M;
    unsigned int m_E;
    unsigned int m_g;
    unsigned int *m_vardegree;
    unsigned int *m_chkdegree;
    unsigned int *m_interleaver;
    double *m_c2v;
    double *m_v2c;
    int *m_hard;
    double *m_Lsyn;

    bool m_clearmsg;
    unsigned int m_method;

    bool m_decloaded;

    unsigned int m_corrveclen;
    double *m_corrvec;

    unsigned int decodeSPA(double *Lch, double *Lapp, unsigned int maxit);
    unsigned int decodeMSA(double *Lch, double *Lapp, unsigned int maxit);

    double L2mutual(double *L, long N);

    void allocateMem();
    void freeMem();
};

#endif

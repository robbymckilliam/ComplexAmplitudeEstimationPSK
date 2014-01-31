#ifndef LDPCDEC_H
#define LDPCDEC_H

#define DECODER_SPA		0
#define DECODER_MSA		1

/// \brief		Main class for LDPC decoding
/// \author		Gottfried Lechner (gottfried.lechner@unisa.edu.au)
/// \version	3.0
/// \date     Created      : March 2008
/// \date     Last modified: April 2012
///
/// This class holds the LDPC decoder.

class CLDPCDec {
public:
    CLDPCDec(const char *filename);
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

    unsigned int getN();
    unsigned int getM();
    unsigned int getK();
    unsigned int getE();

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

///** Wraps a CLDPC just for encoding */
//class LDPCEncoder {
//    public:
//    LDPCEncoder(const char* filename) : 
//        codec(filename),
//        codeword((unsigned int*) malloc(codec->getN() * sizeof (unsigned int))),
//        info(codeword),
//        parity(codeword+codec->getK()),
//        N(codec->getN()),
//        K(codec->getK()),
//        bits(codec->getN())
//    {  
//    }
//        
//    ~LDPCEncoder() {
//        free(codeword);
//    }
//    
//    const std::vector<unsigned int>& encode(const std::vector<unsigned int>& b) {
//        if(b.size() != K) throw "Wrong number of transmit bits";
//        for(int i = 0; i < K; i++) info[i] = b[i];
//        codec->encodeRA(info, parity);
//        for(int i = 0; i < N; i++) bits[i] = codeword[i];
//    }
//    
//    protected:
//     const int N, K;
//     std::vector<unsigned int> bits;
//     std::vector<unsigned int> codeword;
//     const CLDPCDec codec;
//     const unsigned int *codeword, *info, *parity;
//};

#endif

#ifndef TGSW_FUNCTIONS_H
#define TGSW_FUNCTIONS_H

///@file
///@brief This file contains the declaration of TGSW related functions

#include "tfhe_core.h"
#include "tgsw.h"
#include "tlwe.h"
#include "lwebootstrappingkey.h"

// Gsw
EXPORT void tGswKeyGen(TGswKey* result);
EXPORT void tGswSymEncrypt(TGswSample* result, const IntPolynomial* message, double alpha, const TGswKey* key);
EXPORT void tGswSymEncryptInt(TGswSample* result, const int message, double alpha, const TGswKey* key);
EXPORT void tGswSymDecrypt(IntPolynomial* result, const TGswSample* sample, const TGswKey* key, const int Msize);
EXPORT int tGswSymDecryptInt(const TGswSample* sample, const TGswKey* key); 
//do we really decrypt Gsw samples?

// support Functions for TGsw
// Result = 0
EXPORT void tGswClear(TGswSample* result, const TGswParams* params);
// Result += H
EXPORT void tGswAddH(TGswSample* result, const TGswParams* params);
// Result += mu*H
EXPORT void tGswAddMuH(TGswSample* result, const IntPolynomial* message, const TGswParams* params);
// Result += mu*H, mu integer
EXPORT void tGswAddMuIntH(TGswSample* result, const int message, const TGswParams* params);
// Result = tGsw(0)
EXPORT void tGswEncryptZero(TGswSample* result, double alpha, const TGswKey* key);

//fonction de decomposition
EXPORT void tGswTLweDecompH(IntPolynomial* result, const TLweSample* sample, const TGswParams* params);

EXPORT void tGswTorus32PolynomialDecompH(IntPolynomial* result, const TorusPolynomial* sample, const TGswParams* params);
EXPORT void tGswTLweDecompH(IntPolynomial* result, const TLweSample* sample,const TGswParams* params);	

//TODO: Ilaria.Theoreme3.5
EXPORT void tGswExternProduct(TLweSample* result, const TGswSample* a, const TLweSample* b, const TGswParams* params);

// result=result+ (X^ai-1)*bi (ligne 5 de l'algo)
EXPORT void tGswMulByXaiMinusOne(TGswSample* result, int ai, const TGswSample* bk, const TGswParams* params);

//ligne 5 algo,mult externe
EXPORT void tGswExternMulToTLwe(TLweSample* accum, const TGswSample* sample,const TGswParams* params);

/** result = (0,mu) */
EXPORT void tGswNoiselessTrivial(TGswSample* result, const TorusPolynomial* mu, const TGswParams* params);

/** result = result + sample */
EXPORT void tGswAddTo(TGswSample* result, const TGswSample* sample, const TGswParams* params);

/** result = result - sample */
//EXPORT void tGswSubTo(TLweSample* result, const TLweSample* sample, const TLweParams* params);
/** result = result + p.sample */
//EXPORT void tGswAddMulTo(TLweSample* result, int p, const TLweSample* sample, const TLweParams* params);
/** result = result - p.sample */
//EXPORT void tGswSubMulTo(TLweSample* result, int p, const TLweSample* sample, const TLweParams* params);

EXPORT void tfhe_createLweBootstrappingKey(LweBootstrappingKey* bk, const LweKey* key_in, const TGswKey* rgsw_key);
EXPORT void tfhe_bootstrap(LweSample* result, const LweBootstrappingKey* bk, Torus32 mu, const LweSample* x);

EXPORT void tGswToFFTConvert(TGswSampleFFT* result, const TGswSample* source, const TGswParams* params);
EXPORT void tGswFromFFTConvert(TGswSample* result, const TGswSampleFFT* source, const TGswParams* params);
EXPORT void tGswFFTAddH(TGswSampleFFT* result, const TGswParams* params);
EXPORT void tGswFFTClear(TGswSampleFFT* result, const TGswParams* params);
EXPORT void tGswLagrangeHalfCPolynomialDecompH(LagrangeHalfCPolynomial* reps, const LagrangeHalfCPolynomial* pol, const TGswParams* params);
EXPORT void tGswFFTExternMulToTLwe(TLweSampleFFT* accum, TGswSampleFFT* gsw, const TGswParams* params);
EXPORT void tGswFFTMulByXaiMinusOne(TGswSampleFFT* result,TGswSampleFFT* result2, TGswSampleFFT* result3, TGswSampleFFT* result4,TGswSampleFFT* result5,TGswSampleFFT* result6, TGswSampleFFT* result7, TGswSampleFFT* result8, const int ai, const int ai2,const int ai3, const LweBootstrappingKeyFFT* bk, const TGswParams* params, int i);

EXPORT void tfhe_createLweBootstrappingKeyFFT(LweBootstrappingKeyFFT* bk, const LweKey* key_in, const TGswKey* rgsw_key);
EXPORT void tfhe_bootstrapFFT(LweSample* result, const LweBootstrappingKeyFFT* bk, Torus32 mu1, Torus32 mu0, const LweSample* x);

#endif //TGSW_FUNCTIONS_H

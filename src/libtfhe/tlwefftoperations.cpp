#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include <ccomplex>
#include "tfhe_core.h"
#include "numeric_functions.h"
#include "lweparams.h"
#include "lwekey.h"
#include "lwesamples.h"
#include "lwe-functions.h"
#include "tlwe_functions.h"
#include "tgsw_functions.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"
#include "lwebootstrappingkey.h"

using namespace std;

#ifndef NDEBUG
extern const TLweKey* debug_accum_key;
extern const LweKey* debug_extract_key;
extern const LweKey* debug_in_key;
#endif


EXPORT void tLweToFFTConvert(TLweSampleFFT* result, const TLweSample* source, const TLweParams* params){
    const int k = params->k;
    
    for (int i = 0; i <= k; ++i)
	TorusPolynomial_ifft(result->a+i,source->a+i);
}

EXPORT void tLweFromFFTConvert(TLweSample* result, const TLweSampleFFT* source, const TLweParams* params){
    const int k = params->k;
    
    for (int i = 0; i <= k; ++i)
	TorusPolynomial_fft(result->a+i,source->a+i);
}



//Arithmetic operations on TLwe samples
/** result = (0,0) */
EXPORT void tLweFFTClear(TLweSampleFFT* result, const TLweParams* params){
    int k = params->k;

    for (int i = 0; i <= k; ++i) LagrangeHalfCPolynomialClear(&result->a[i]);
    result->current_variance = 0.;
}

/** result = (0,mu) */
EXPORT void tLweFFTNoiselessTrivial(TLweSampleFFT* result, const TorusPolynomial* mu, const TLweParams* params){
    int k = params->k;

    for (int i = 0; i < k; ++i) LagrangeHalfCPolynomialClear(&result->a[i]);
    TorusPolynomial_ifft(result->b, mu);
    result->current_variance = 0.;
}

/** result = (0,mu) where mu is constant*/
EXPORT void tLweFFTNoiselessTrivialT(TLweSampleFFT* result, const Torus32 mu, const TLweParams* params){
    const int k = params->k;
    
    for (int i = 0; i < k; ++i) 
	LagrangeHalfCPolynomialClear(&result->a[i]);
    LagrangeHalfCPolynomialSetTorusConstant(result->b,mu);
    result->current_variance = 0.;
}

/** result = result + sample */
EXPORT void tLweFFTAddTo(TLweSampleFFT* result, const TLweSampleFFT* sample, const TLweParams* params);
//Let's postpone this to make sure we actually need it
//{
//    int k = params->k;
//
//    for (int i = 0; i < k; ++i) 
//	AddToLagrangeHalfCPolynomial(&result->a[i], &sample->a[i]);
//    AddToLagrangeHalfCPolynomial(result->b, sample->b);
//    result->current_variance += sample->current_variance; //à revoir//OK si c'est la variance
//}

/** result = result - sample */
EXPORT void tLweFFTSubTo(TLweSampleFFT* result, const TLweSampleFFT* sample, const TLweParams* params);

/** result = result + p.sample */
EXPORT void tLweFFTAddMulZTo(TLweSampleFFT* result, int p, const TLweSampleFFT* sample, const TLweParams* params);
//Let's postpone this to make sure we actually need it
//{
//    int k = params->k;
//
//    for (int i = 0; i < k; ++i) 
//	torusPolynomialAddMulZTo(&result->a[i], p, &sample->a[i]);
//    torusPolynomialAddMulZTo(result->b, p, sample->b);
//    result->current_variance += (p*p)*sample->current_variance;
//}

/** result = result - p.sample */
EXPORT void tLweFFTSubMulZTo(TLweSampleFFT* result, int p, const TLweSampleFFT* sample, const TLweParams* params);


EXPORT void tLweFFTAddMulRTo(TLweSampleFFT* result, const LagrangeHalfCPolynomial* p, const TLweSampleFFT* sample, const TLweParams* params) {
    const int k = params->k;
    
    for (int i=0; i<=k; i++)
	LagrangeHalfCPolynomialAddMul(result->a+i,p,sample->a+i);
}

EXPORT void tLweFFTMulR(TLweSampleFFT* result, const LagrangeHalfCPolynomial* p, const TLweSampleFFT* sample, const TLweParams* params) {
    const int k = params->k;
    
    for (int i=0; i<=k; i++)
	LagrangeHalfCPolynomialMul(result->a+i,p,sample->a+i);
}

EXPORT void tLweFFTSubMulRTo(TLweSampleFFT* result, const LagrangeHalfCPolynomial* p, const TLweSampleFFT* sample, const TLweParams* params) {
    const int k = params->k;
    
    for (int i=0; i<=k; i++)
	LagrangeHalfCPolynomialSubMul(result->a+i,p,sample->a+i);
}

    
EXPORT void tGswToFFTConvert(TGswSampleFFT* result, const TGswSample* source, const TGswParams* params) {
    const int kpl = params->kpl;
    
    for (int p=0; p<kpl; p++)
	tLweToFFTConvert(result->all_samples+p, source->all_sample+p, params->tlwe_params);
}

EXPORT void tGswFromFFTConvert(TGswSample* result, const TGswSampleFFT* source, const TGswParams* params){
    const int kpl = params->kpl;
    
    for (int p=0; p<kpl; p++)
	tLweFromFFTConvert(result->all_sample+p, source->all_samples+p, params->tlwe_params);
}

EXPORT void tGswFFTAddH(TGswSampleFFT* result, const TGswParams* params) {
    const int k = params->tlwe_params->k;
    const int l = params->l;

    for (int j=0; j<l; j++) {
    	Torus32 hj = params->h[j];
    	for (int i=0; i<=k; i++)
	   LagrangeHalfCPolynomialAddTorusConstant(&result->sample[i][j].a[i],hj); 
    }

}

/*EXPORT void tGswFFTAddH2(TGswSampleFFT* result,TGswSampleFFT* result2, const TGswParams* params) {
    const int k = params->tlwe_params->k;
    const int l = params->l;

    for (int j=0; j<l; j++) {

    	for (int i=0; i<=k; i++)
	   LagrangeHalfCPolynomialAddTorusConstant2(&result->sample[i][j].a[i],&result2->sample[i][j].a[i]); 
    }

}
*/
EXPORT void tGswFFTAddH8(TGswSampleFFT* result,TGswSampleFFT* result2,TGswSampleFFT* result3,TGswSampleFFT* result4,TGswSampleFFT* result5,TGswSampleFFT* result6,TGswSampleFFT* result7,TGswSampleFFT* result8, const TGswParams* params) {
    const int k = params->tlwe_params->k;
    const int l = params->l;

    for (int j=0; j<l; j++) {

    	for (int i=0; i<=k; i++)
	   LagrangeHalfCPolynomialAddTorusConstant8(&result->sample[i][j].a[i],&result2->sample[i][j].a[i],&result3->sample[i][j].a[i],&result4->sample[i][j].a[i],&result5->sample[i][j].a[i],&result6->sample[i][j].a[i],&result7->sample[i][j].a[i],&result8->sample[i][j].a[i]); 
    }

}

EXPORT void tGswFFTClear(TGswSampleFFT* result, const TGswParams* params) {
    const int kpl = params->kpl;

    for (int p=0; p<kpl; p++)
	tLweFFTClear(result->all_samples+p, params->tlwe_params);
}    

EXPORT void tGswLagrangeHalfCPolynomialDecompH(LagrangeHalfCPolynomial* reps, const LagrangeHalfCPolynomial* pol, const TGswParams* params) {
    const int l = params->l;
    const int N = params->tlwe_params->N;
    //TODO attention, this prevents parallelization...
    static TorusPolynomial* a = new_TorusPolynomial(N);
    static IntPolynomial* deca = new_IntPolynomial_array(l,N);

    TorusPolynomial_fft(a,pol);
    tGswTorus32PolynomialDecompH(deca, a, params);
    for (int j=0; j<l; j++) {
	IntPolynomial_ifft(reps+j,deca+j);
    }
}

EXPORT void tGswFFTExternMulToTLwe(TLweSampleFFT* accum, TGswSampleFFT* gsw, const TGswParams* params) {
    const TLweParams* tlwe_params=params->tlwe_params;
    const int k = tlwe_params->k;
    const int l = params->l;
    const int kpl = params->kpl;
    const int N = tlwe_params->N;
    //TODO attention, this prevents parallelization...
    static LagrangeHalfCPolynomial* decomps=new_LagrangeHalfCPolynomial_array(kpl,N);

    for (int i=0; i<=k; i++)
	tGswLagrangeHalfCPolynomialDecompH(decomps+i*l,accum->a+i, params);
    tLweFFTClear(accum, tlwe_params);
    for (int p=0; p<kpl; p++)
	tLweFFTAddMulRTo(accum, decomps+p, gsw->all_samples+p, tlwe_params);
}

//修改
EXPORT void tGswFFTMulByXaiMinusOne(TGswSampleFFT* result,TGswSampleFFT* result2, TGswSampleFFT* result3, TGswSampleFFT* result4,TGswSampleFFT* result5,TGswSampleFFT* result6, TGswSampleFFT* result7, TGswSampleFFT* result8, const int ai, const int ai2,const int ai3, const LweBootstrappingKeyFFT* bk, const TGswParams* params, int i) {
    const TLweParams* tlwe_params=params->tlwe_params;
    const int k = tlwe_params->k;
    //const int l = params->l;
    const int kpl = params->kpl;
    const int N = tlwe_params->N;
    //on calcule x^ai-1 en fft
    //TODO attention, this prevents parallelization...
    static LagrangeHalfCPolynomial* xaim1=new_LagrangeHalfCPolynomial(N);
    static LagrangeHalfCPolynomial* xaim2=new_LagrangeHalfCPolynomial(N);
    static LagrangeHalfCPolynomial* xaim3=new_LagrangeHalfCPolynomial(N);
    static LagrangeHalfCPolynomial* xaim4=new_LagrangeHalfCPolynomial(N);
    static LagrangeHalfCPolynomial* xaim5=new_LagrangeHalfCPolynomial(N);
    static LagrangeHalfCPolynomial* xaim6=new_LagrangeHalfCPolynomial(N);
    static LagrangeHalfCPolynomial* xaim7=new_LagrangeHalfCPolynomial(N);

    //--lagrangehalfc_impl.cpp(fft--fftw)--
    //--lagrangehalfc_arithmetic(include)
/*
    LagrangeHalfCPolynomialSetXaiMinusOne(xaim1,ai);
    LagrangeHalfCPolynomialSetXaiMinusOne(xaim2,ai2);
    LagrangeHalfCPolynomialSetXaiMinusOne2(xaim3,ai,ai2);
*/


    LagrangeHalfCPolynomialSetXaiMinusOne3(xaim1,ai,ai2,ai3);
    LagrangeHalfCPolynomialSetXaiMinusOne2(xaim2,ai,ai2);
    LagrangeHalfCPolynomialSetXaiMinusOne2(xaim3,ai,ai3);
    LagrangeHalfCPolynomialSetXaiMinusOne2(xaim4,ai2,ai3);
    LagrangeHalfCPolynomialSetXaiMinusOne(xaim5,ai);
    LagrangeHalfCPolynomialSetXaiMinusOne(xaim6,ai2);
    LagrangeHalfCPolynomialSetXaiMinusOne(xaim7,ai3);



    for (int p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bk->bk[i].all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result->all_samples[p].a;

        const LagrangeHalfCPolynomial* in_s2 =bk->bk[i+1].all_samples[p].a;
        LagrangeHalfCPolynomial* out_s2 = result2->all_samples[p].a;

        const LagrangeHalfCPolynomial* in_s3 =bk->bk[i+2].all_samples[p].a;
        LagrangeHalfCPolynomial* out_s3 = result3->all_samples[p].a;

        const LagrangeHalfCPolynomial* in_s4 =bk->bk[i+3].all_samples[p].a;
        LagrangeHalfCPolynomial* out_s4 = result4->all_samples[p].a;

        const LagrangeHalfCPolynomial* in_s5 =bk->bk[i+4].all_samples[p].a;
        LagrangeHalfCPolynomial* out_s5 = result5->all_samples[p].a;

        const LagrangeHalfCPolynomial* in_s6 =bk->bk[i+5].all_samples[p].a;
        LagrangeHalfCPolynomial* out_s6 = result6->all_samples[p].a;

        const LagrangeHalfCPolynomial* in_s7 =bk->bk[i+6].all_samples[p].a;
        LagrangeHalfCPolynomial* out_s7 = result7->all_samples[p].a;


        for (int j=0; j<=k; j++)
{
            LagrangeHalfCPolynomialMul(&out_s[j], xaim1, &in_s[j]); 
            LagrangeHalfCPolynomialMul(&out_s2[j], xaim2, &in_s2[j]); 
            LagrangeHalfCPolynomialMul(&out_s3[j], xaim3, &in_s3[j]); 
            LagrangeHalfCPolynomialMul(&out_s4[j], xaim4, &in_s4[j]); 
            LagrangeHalfCPolynomialMul(&out_s5[j], xaim5, &in_s5[j]);
            LagrangeHalfCPolynomialMul(&out_s6[j], xaim6, &in_s6[j]);
            LagrangeHalfCPolynomialMul(&out_s7[j], xaim7, &in_s7[j]); 
}
    }

/*    for (int p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bki2->all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result2->all_samples[p].a;


        for (int j=0; j<=k; j++)
            LagrangeHalfCPolynomialMul(&out_s[j], xaim2, &in_s[j]); 
    }

    for (int p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bki3->all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result3->all_samples[p].a;
        for (int j=0; j<=k; j++)
            LagrangeHalfCPolynomialMul(&out_s[j], xaim3, &in_s[j]); 
    }

    for (int p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bki4->all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result4->all_samples[p].a;
        for (int j=0; j<=k; j++)
            LagrangeHalfCPolynomialMul(&out_s[j], xaim4, &in_s[j]); 
    }

    for (int p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bki5->all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result5->all_samples[p].a;
        for (int j=0; j<=k; j++)
            LagrangeHalfCPolynomialMul(&out_s[j], xaim5, &in_s[j]); 
    }

    for (int p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bki6->all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result6->all_samples[p].a;
        for (int j=0; j<=k; j++)
            LagrangeHalfCPolynomialMul(&out_s[j], xaim6, &in_s[j]); 
    }

    for (int p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bki7->all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result7->all_samples[p].a;
        for (int j=0; j<=k; j++)
            LagrangeHalfCPolynomialMul(&out_s[j], xaim7, &in_s[j]); 
    }
*/
}

EXPORT void tfhe_createLweBootstrappingKeyFFT(
	LweBootstrappingKeyFFT* bk, 
	const LweKey* key_in, 
	const TGswKey* rgsw_key) {
    assert(bk->bk_params==rgsw_key->params);
    assert(bk->in_out_params==key_in->params);

    const LweParams* in_out_params = bk->in_out_params; 
    const TGswParams* bk_params = bk->bk_params;
    const TLweParams* accum_params = bk_params->tlwe_params;
    const LweParams* extract_params = &accum_params->extracted_lweparams;
    
    //LweKeySwitchKey* ks; ///< the keyswitch key (s'->s)
    const TLweKey* accum_key = &rgsw_key->tlwe_key;
    LweKey* extracted_key = new_LweKey(extract_params);
    tLweExtractKey(extracted_key, accum_key);
    lweCreateKeySwitchKey(bk->ks, extracted_key, key_in);
    delete_LweKey(extracted_key);
    
    //TGswSample* bk; ///< the bootstrapping key (s->s")
    TGswSample* tmpsample = new_TGswSample(bk_params);
    int* kin = key_in->key;
    const double alpha = accum_params->alpha_min;
    const int n = in_out_params->n;
    //cout << "tfhe_createLweBootstrappingKeyFFT::tGswSymEncryptInt _LLF" << endl;
//    for (int i=0; i<(n/3)*3+3; i=i+3) {
//    for(int j=(n/3)*3; j<n+3; j=j+1)
//        {kin[j]=0;}

    for (int i=0; i<n; i=i+3) {

	
	//cout << "starting the test..." << i <<kin[i]<< endl;
	//bk1
	tGswSymEncryptInt(tmpsample, kin[i]*kin[i+1]*kin[i+2], alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3], tmpsample, bk_params);

	tGswSymEncryptInt(tmpsample, (1-kin[i]*kin[i+1]*(kin[i+2]-1)), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3+1], tmpsample, bk_params);

	//bk2
	tGswSymEncryptInt(tmpsample, (1-kin[i]*(kin[i+1]-1)*kin[i+2]), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3+2], tmpsample, bk_params);

	tGswSymEncryptInt(tmpsample, (1-(kin[i]-1)*kin[i+1]*kin[i+2]), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3+3], tmpsample, bk_params);
	
	//bk3
	tGswSymEncryptInt(tmpsample, (kin[i]*(kin[i+1]-1)*(kin[i+2]-1)), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3+4], tmpsample, bk_params);

	tGswSymEncryptInt(tmpsample, ((kin[i]-1)*kin[i+1]*(kin[i+2]-1)), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3+5], tmpsample, bk_params);

	//bk4
	tGswSymEncryptInt(tmpsample, ((kin[i]-1)*(kin[i+1]-1)*kin[i+2]), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3+6], tmpsample, bk_params);

	tGswSymEncryptInt(tmpsample, (1-(kin[i]-1)*(kin[i+1]-1)*(kin[i+2]-1)), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i*8/3+7], tmpsample, bk_params);
    }

/*
    //这里设置i不整除3，然后对最后几个处理
	const int i=(n/3)*3;
	cout << "i..." << i << endl;
	tGswSymEncryptInt(tmpsample, kin[i]*kin[i+1], alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i], tmpsample, bk_params);

	tGswSymEncryptInt(tmpsample, (1-(kin[i]-1)*kin[i+1]), alpha, rgsw_key);
	tGswToFFTConvert(&bk->bk[i+1], tmpsample, bk_params);

	tGswSymEncryptInt(tmpsample, (1-(kin[i+1]-1)*kin[i]), alpha, rgsw_key);
	tGswToFFTConvert(&bk2->bk[i], tmpsample, bk_params);

	tGswSymEncryptInt(tmpsample, (kin[i+1]-1)*(kin[i]-1), alpha, rgsw_key);
	tGswToFFTConvert(&bk2->bk[i+1], tmpsample, bk_params);		
*/
    delete_TGswSample(tmpsample);
}


EXPORT void tfhe_bootstrapFFT(LweSample* result, const LweBootstrappingKeyFFT* bk,  Torus32 mu1, Torus32 mu0, const LweSample* x){
    const Torus32 ab=(mu1+mu0)/2;
    const Torus32 aa = mu0-ab; // aa=(mu1-mu0)/2;
    const TGswParams* bk_params = bk->bk_params;
    const TLweParams* accum_params = bk_params->tlwe_params;
    const LweParams* extract_params = &accum_params->extracted_lweparams;
    const LweParams* in_out_params = bk->in_out_params;
    const int n=in_out_params->n;
    const int N=accum_params->N;
    const int Ns2=N/2;
    const int Nx2= 2*N;
    
    
    TorusPolynomial* testvect=new_TorusPolynomial(N);
    TorusPolynomial* testvectbis=new_TorusPolynomial(N);


    int barb=modSwitchFromTorus32(x->b,Nx2);
    //je definis le test vector (multiplié par a inclus !
    for (int i=0;i<Ns2;i++)
       testvect->coefsT[i]=aa;
    for (int i=Ns2;i<N;i++)
       testvect->coefsT[i]=-aa;
    torusPolynomialMulByXai(testvectbis, barb, testvect);

    // Accumulateur 
    TLweSample* acc = new_TLweSample(accum_params);
    TLweSampleFFT* accFFT = new_TLweSampleFFT(accum_params);

    // acc and accFFt will be used for tfhe_bootstrapFFT, acc1=acc will be used for tfhe_bootstrap
    tLweNoiselessTrivial(acc, testvectbis, accum_params);
    tLweToFFTConvert(accFFT, acc, accum_params);

    TGswSample* temp = new_TGswSample(bk_params);

    //修改
    TGswSampleFFT* tempFFT = new_TGswSampleFFT(bk_params);
    TGswSampleFFT* tempFFT2 = new_TGswSampleFFT(bk_params);
    TGswSampleFFT* tempFFT3 = new_TGswSampleFFT(bk_params);
    TGswSampleFFT* tempFFT4 = new_TGswSampleFFT(bk_params);
    TGswSampleFFT* tempFFT5 = new_TGswSampleFFT(bk_params);
    TGswSampleFFT* tempFFT6 = new_TGswSampleFFT(bk_params);
    TGswSampleFFT* tempFFT7 = new_TGswSampleFFT(bk_params);
    TGswSampleFFT* tempFFT8 = new_TGswSampleFFT(bk_params);

//NICOLAS: j'ai ajouté ce bloc
#ifndef NDEBUG
    TorusPolynomial* phase = new_TorusPolynomial(N);
    int correctOffset = barb;
    cout << "starting the test..." << endl;
#endif
    // the index 1 is given when we don't use the fft
    for (int i=0; i<n; i=i+3) {
//    	clock_t begin = clock();
        int bara=modSwitchFromTorus32(-x->a[i],Nx2);
	int bara2=modSwitchFromTorus32(-x->a[i+1],Nx2);
  	int bara3=modSwitchFromTorus32(-x->a[i+2],Nx2);      
	//--tgswfunciton.h--tlwefftoperations.cppbk[i*8/3+1]
        tGswFFTMulByXaiMinusOne(tempFFT, tempFFT2, tempFFT3, tempFFT4, tempFFT5, tempFFT6, tempFFT7, tempFFT8, bara, bara2, bara3, bk, bk_params,i*8/3);

/*
	tGswFFTAddH2(tempFFT, tempFFT2, bk_params);
	tGswFFTAddH2(tempFFT, tempFFT3, bk_params);
	tGswFFTAddH2(tempFFT, tempFFT4, bk_params);
	tGswFFTAddH2(tempFFT, tempFFT5, bk_params);
	tGswFFTAddH2(tempFFT, tempFFT6, bk_params);
	tGswFFTAddH2(tempFFT, tempFFT7, bk_params);

	tGswFFTAddH2(tempFFT, bk4->bk+i+1, bk_params);

*/
	tGswFFTAddH8(tempFFT, tempFFT2, tempFFT3, tempFFT4, tempFFT5, tempFFT6, tempFFT7, &bk->bk[i*8/3+7],bk_params);
 //   clock_t end = clock();
 //   cout << "finished bootstrapping (microsecs)1... " << (end-begin) << endl;


        tGswFFTExternMulToTLwe(accFFT, tempFFT, bk_params);
//    clock_t end2 = clock();
//    cout << "finished bootstrapping (microsecs)2... " << (end2-end) << endl;

//NICOLAS: et surtout, j'ai ajouté celui-ci!
#ifndef NDEBUG
            tLweFromFFTConvert(acc, accFFT, accum_params);
        tLwePhase(phase,acc,debug_accum_key);  //celui-ci, c'est la phase de acc (FFT)
	if (debug_in_key->key[i]==1) correctOffset = (correctOffset+bara)%Nx2; 
        torusPolynomialMulByXai(testvectbis, correctOffset, testvect); //celui-ci, c'est la phase idéale (calculée sans bruit avec la clé privée)
	for (int j=0; j<N; j++) {
	       printf("Iteration %d, index %d: phase %d vs noiseless %d\n",i,j,phase->coefsT[j], testvectbis->coefsT[j]);
	}
#endif

    }
    tLweFromFFTConvert(acc, accFFT, accum_params);


    LweSample* u = new_LweSample(extract_params);
    tLweExtractLweSample(u, acc, extract_params, accum_params);
    u->b += ab;
    
    lweKeySwitch(result, bk->ks, u);
    


    delete_LweSample(u);
    delete_TGswSampleFFT(tempFFT); 
    delete_TGswSample(temp);
    delete_TLweSampleFFT(accFFT);
    delete_TLweSample(acc);
    delete_TorusPolynomial(testvectbis);
    delete_TorusPolynomial(testvect);
}

#ifndef LAGRANGEHALFC_ARITHMETIC_H
#define LAGRANGEHALFC_ARITHMETIC_H

///@file
///@brief This file contains the declaration of operations on LagrangeHalfC polynomials

#include "polynomials.h"

//initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ constructor)
EXPORT void init_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj, const int N);

//destroys the LagrangeHalfCPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj);


/**
 * FFT functions 
 */
EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p);
EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p);
EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p);

//MISC OPERATIONS
/** sets to zero */
EXPORT void LagrangeHalfCPolynomialClear(LagrangeHalfCPolynomial* result);

/** sets to this torus32 constant */
EXPORT void LagrangeHalfCPolynomialSetTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu);
EXPORT void LagrangeHalfCPolynomialAddTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 cst);
EXPORT void LagrangeHalfCPolynomialAddTorusConstant2(LagrangeHalfCPolynomial* result, LagrangeHalfCPolynomial* result2);
EXPORT void LagrangeHalfCPolynomialAddTorusConstant8(LagrangeHalfCPolynomial* result, LagrangeHalfCPolynomial* result2,LagrangeHalfCPolynomial* result3, LagrangeHalfCPolynomial* result4,LagrangeHalfCPolynomial* result5, LagrangeHalfCPolynomial* result6,LagrangeHalfCPolynomial* result7, LagrangeHalfCPolynomial* result8);
/** sets to X^ai-1 */
EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne2(LagrangeHalfCPolynomial* result, const int ai, const int ai2);
EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne(LagrangeHalfCPolynomial* result, const int ai);
EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne3(LagrangeHalfCPolynomial* result, const int ai, const int ai2, const int ai3);
/** multiplication via direct FFT */
EXPORT void torusPolynomialMultFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);
EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);
EXPORT void torusPolynomialSubMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2);

/** termwise multiplication in Lagrange space */
EXPORT void LagrangeHalfCPolynomialMul(
	LagrangeHalfCPolynomial* result, 
	const LagrangeHalfCPolynomial* a, 
	const LagrangeHalfCPolynomial* b);

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialClear(
	LagrangeHalfCPolynomial* reps);

EXPORT void LagrangeHalfCPolynomialAddTo(
	LagrangeHalfCPolynomial* accum, 
	const LagrangeHalfCPolynomial* a);

EXPORT void LagrangeHalfCPolynomialAddMul(
	LagrangeHalfCPolynomial* accum, 
	const LagrangeHalfCPolynomial* a, 
	const LagrangeHalfCPolynomial* b);

EXPORT void LagrangeHalfCPolynomialSubMul(
	LagrangeHalfCPolynomial* accum, 
	const LagrangeHalfCPolynomial* a, 
	const LagrangeHalfCPolynomial* b);

#endif // LAGRANGEHALFC_ARITHMETIC_H

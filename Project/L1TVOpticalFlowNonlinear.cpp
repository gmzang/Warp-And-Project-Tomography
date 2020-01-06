//
//
//#include"L1TVOpticalFlow.h"
//#include<iostream>
//#include <algorithm>
//#include <string.h>
////#include "math.h"
//
//#include <omp.h>
//#include "tools.h"
////#include "linearInterpolation.h"
////#include "cubicInterpolation.h"
//#include <stdio.h>
//#include <sys/types.h>
//
//#include<string_utils.h>
//#include "advection.h"
//
////#include"linear_solver.h"
//
//#include "Vector3.h"
//
//
////#include "advection2D.h"
//
//#include <cstddef>
////#include "xmalloc.c"
//#include <ctime>
//# include <omp.h>
//
//#define DIMENSION_SIZE 6
//int ROI[] = { 10,8,5,50,5,5 };
//
////#include "vec3.h"
////using namespace geom;
////typedef geom::vec3<float>	real3;
////#include "CImg.h"
////using namespace cimg_library;
////using namespace std;
////namespace cl = cimg_library;
////
////typedef float				real;
//////typedef geom::vec3<real>	real3;
////typedef cl::CImg<real>		image;
////typedef Vector3	real3;
//
//int interpo = 5;
//using namespace std;
////real	m_aabb[6];// = { -1000,1000,-1000,1000,-1200,1200 };
//
//template< typename real, typename image >
//image blur_and_downsample(const image &I, const int nx, const int ny, const int nz, const real sigma) {
//	return I.get_blur(sigma).get_resize(nx, ny, nz, I.spectrum(), 5);
//}
//
////template< typename real, typename image >
////image blur_and_downsample(const image &I, const real eta, const real sigma) {
////	return blur_and_downsample(I, I.width()*eta, I.height()*eta, I.depth()*eta, sigma);
////}
//template< typename real, typename image >
//image blur_and_downsample(const image &I, const real eta, const real sigma) {
//	return blur_and_downsample(I, I.width()*eta, I.height()*eta, I.depth(), sigma);
//}
//template< typename real, typename image >
//image Y_blur_and_downsample(const image &I, const real eta, const real sigma) {
//	return blur_and_downsample(I, I.width()*eta, I.height()*eta, I.depth(), sigma);
//}
//
//
////void doWarp(const float *image1, const float *image2, float *v1, float *v2, float *ut, float *ux, float *uy, float *uxt, float *uyt, float *uxx, float *uxy, float *uyx, float *uyy, const int *sizeMat);
////void doWarp( float *image1,  float *image2, float *v1, float *v2, float *ut, float *ux, float *uy,   int *sizeMat);
////void doWarp3D( float *image1,  float *image2, float *v1, float *v2, float *v3, float *ut, float *ux, float *uy, float *uz,  int *sizeMat);
//
////void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//// [x(:,:,j,1),x(:,:,j,2),yRes] = ((uTmp(:,:,j),uTmp(:,:,j+1),tol,alpha,maxIt,typeNorm,
////     xT,yT,stepsize,discretization,numberOfWarps,huberEpsilon,gradientConstancy)']);
//
//
////optical_flow_horn_schunk_multiscale
//#if DIMENSION_SIZE ==2
//
//void L1TVOpticalFlowNonlinear_MultiScale(const int *dim, const image &image1, const image &image2,  image &_inputV,
//										 image &_inputY1, image &_inputY2, image &_inputY3, float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps,
//										float eta, float sigma, int scales)
//{
//	//if (levels > 1)
//	//{
//
//	//	std::cout << "V size: " << _inputV.width() << " " << _inputV.height() << std::endl;
//	//	image I1_coarse = blur_and_downsample(image1, eta, sigma);
//	//	image I2_coarse = blur_and_downsample(image2, eta, sigma);
//	//	image V_coarse = Y_blur_and_downsample(_inputV, eta, sigma);
//	//	// maybe need to seperate the y1, y2, y
//	//	//image Y_coarse= blur_and_downsample(_inputY, eta, sigma);
//	//	image Y_coarse = Y_blur_and_downsample(_inputY, eta, sigma);
//	//	image V_coarse_init = V_coarse;
//	//	int tdim[] = { I1_coarse.width(), I1_coarse.height() };
//	//	std::cout << "tdim size: " << tdim[0] <<" "<<tdim[1]<< std::endl;
//
//	//	//int tdim[] = { I1_coarse.width(), I1_coarse.height(), I1_coarse.depth() };
//	//	//optical_flow_horn_schunk_multiscale(tdim, aabb, I1_coarse, I2_coarse, uvw_coarse, opts, point_data_coarse, cell_data_coarse, level - 1);
//	//
//
//	//	 L1TVOpticalFlowNonlinear_MultiScale(tdim, I1_coarse, I2_coarse, V_coarse,
//	//										 Y_coarse, _tol, _lambda, _maxIterations, _norm, _numberOfWarps,eta,  sigma,  levels-1);
//
//
//
//	//	 V_coarse += (V_coarse - V_coarse_init).get_resize(dim[0], dim[1], 2,3, 3);
//
//	//}
//	//std::cout <<  "optical flow level: " << levels << std::endl;
//
//
//	//return L1TVOpticalFlowNonlinear(dim, image1, image2, _inputV, _inputY, _tol, _lambda, _maxIterations, _norm,
//	//								_numberOfWarps);
//	
//	
//		//image I1_coarse = blur_and_downsample(image1, eta, sigma);
//		//image I2_coarse = blur_and_downsample(image2, eta, sigma);
//		//image V_coarse = Y_blur_and_downsample(_inputV, eta, sigma);
//		image I1_blur = image1.get_blur(sigma);
//		image I2_blur = image2.get_blur(sigma);
//
//		image tempI1 = I1_blur;
//		image tempI2 = I2_blur;
//		image tempV= _inputV;
//		image tempY1 = _inputY1;
//		image tempY2 = _inputY2;
//		image tempY3 = _inputY3;
//
//		//image I3_blur = image3.get_blur(sigma);
//		//float **Im1s;
//		//float **Im2s;
//		int size = image1.width()*image1.height();
//		//Im1s[0] = new float[size];
//		//Im2s[0] = new float[size];
//		imagelist Im1s;
//		imagelist Im2s;
//		imagelist inputVs;
//		imagelist inputY1s;
//		imagelist inputY2s;
//		imagelist inputY3s;
//
//		imagelist dims;
//		int twidth = I1_blur.width();
//		int theight = I1_blur.height();
//		int tdepth = I1_blur.depth();
//		float t_eta = 1.0f;
//		for (int i = 0; i < scales; i++)
//		{
//
//
//
//			Im1s.insert(I1_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth, 1, 5));
//			Im2s.insert(I2_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth, 1, 5));
//			inputVs.insert(_inputV.get_resize(twidth*t_eta, theight*t_eta, tdepth, 2, 5));
//			inputY1s.insert(_inputY1.get_resize(twidth*t_eta, theight*t_eta, tdepth, 2, 5));
//			inputY2s.insert(_inputY2.get_resize(twidth*t_eta, theight*t_eta, tdepth, 2, 5));
//			inputY3s.insert(_inputY3.get_resize(twidth*t_eta, theight*t_eta, tdepth, 1, 5));
//
//			t_eta = t_eta*eta;
//
//
//
//
//
//	
//		/*	inputY4s.insert(_inputY4.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//			inputY4s.insert(tempY4);
//			image timg7 = tempY4.get_resize(tempY4.width()*eta, tempY4.height()*eta, tempY4.depth()*eta, 1, 5);
//			tempY3 = timg6;*/
//
//			/*Im1s.insert(tempI1);
//				image timg = tempI1.get_resize(tempI1.width()*eta, tempI1.height()*eta, tempI1.depth(), 1, 5);
//				tempI1 = timg;
//
//				Im2s.insert(tempI2);
//				image timg2 = tempI2.get_resize(tempI2.width()*eta, tempI2.height()*eta, tempI2.depth(), 1, 5);
//				tempI2 = timg2;
//
//				inputVs.insert(tempV);
//				image timg3 = tempV.get_resize(tempV.width()*eta, tempV.height()*eta, tempV.depth(), 2, 5);
//				tempV = timg3;
//
//				inputY1s.insert(tempY1);
//				image timg4 = tempY1.get_resize(tempY1.width()*eta, tempY1.height()*eta, tempY1.depth(), 2, 5);
//				tempY1 = timg4;
//
//				inputY2s.insert(tempY2);
//				image timg5 = tempY2.get_resize(tempY2.width()*eta, tempY2.height()*eta, tempY2.depth(), 2, 5);
//				tempY2 = timg5;
//
//				inputY3s.insert(tempY3);
//				image timg6 = tempY3.get_resize(tempY3.width()*eta, tempY3.height()*eta, tempY3.depth(), 1, 5);
//				tempY3 = timg6;
//*/
//		}
//		Im1s.reverse();
//		Im2s.reverse();
//		inputVs.reverse();
//		//inputYs.reverse();
//		inputY1s.reverse();
//		inputY2s.reverse();
//		inputY3s.reverse();
//		cout << "Test 1" << endl;
//		////for (CImgList<>::iterator it = list.begin(); it<list.end(); ++it) (*it).mirror('x');
//		for (int i = 0; i < inputVs.size(); i++)
//		{
//			cout << "At levels " << inputVs.size()-i << endl;
//			int _dim[] = {Im1s(i).width(),Im1s(i).height(),Im1s(i).depth() };
//			int t_dim[] = { inputY1s(i).width(),inputY1s(i).height(),inputY1s(i).depth() };
//			//cout << "_dim[]  :" << _dim[0] << " " << _dim[1] << " " << _dim[2] << endl;
//			//cout << "_dims[]  :" << _dim[0] << " " << _dim[1] << " " << _dim[2] << endl;
//			//L1TVOpticalFlowNonlinear(_dim, Im1s(i), Im2s(i), inputVs(i), inputYs(i), _tol, _lambda, _maxIterations, _norm,
//			//									_numberOfWarps);
//			L1TVOpticalFlowNonlinear(_dim, Im1s(i), Im2s(i), inputVs(i), inputY1s(i), inputY2s(i), inputY3s(i), _tol, _lambda, _maxIterations, _norm,
//				_numberOfWarps);
//
//			cout << "Test 2" << endl;
//
//			if (i < inputVs.size() - 1)
//			{
//				int _dims[] = { inputVs(i + 1).width(), inputVs(i + 1).height(), inputVs(i + 1).depth() };
//
//				inputVs(i+1)= inputVs(i).get_resize(_dims[0], _dims[1], _dims[2], 2, 5);
//				//inputVs(i + 1) = inputVs(i).get_resize(inputVs(i + 1).width(), inputVs(i + 1).height(), inputVs(i + 1).depth(), 2, 5);
//				//inputYs(i + 1) = inputYs(i).get_resize(inputYs(i + 1).width(), inputYs(i + 1).height(), inputYs(i + 1).depth(), 1, 5);
//				inputY1s(i + 1) = inputY1s(i).get_resize(_dims[0], _dims[1], _dims[2], 2, 5);
//				inputY2s(i + 1) = inputY2s(i).get_resize(_dims[0], _dims[1], _dims[2], 2, 5);
//				inputY3s(i + 1) = inputY3s(i).get_resize(_dims[0], _dims[1], _dims[2], 1, 5);
//
//			}
//
//
//			//cout << inputVs.size()<<" "<< inputVs(i).width() << " " << inputVs(i).height() << endl;
//
//		}
//		////cout << Im1s.size() << endl;
//		////system("pause");
//
//}
//
//
//
//
//
//void L1TVOpticalFlowNonlinear( int *dim, const image & image1, const image & image2, image & _inputV, image &_inputY1, image &_inputY2, image &_inputY3,
//	float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps)
//
//{
//
//	//const int maxIterations = (int)mxGetScalar(prhs[4]);
//
//	//const size_t *sizeImage = mxGetDimensions(prhs[0]);
//
//	//float *u1 = image1;
//	//float *u2 = image2;
//	//image u1 = image1;
//	//image u2 = image2;
//
//	const float tol = _tol;
//	const float lambda = _lambda;
//
//	//mexPrintf("tol: %f, lambda: %f; ", tol, lambda);
//
//	const int maxIterations = _maxIterations;
//
//	//const size_t *sizeImage = mxGetDimensions(prhs[0]);
//	//int *sizeImage = new int[2];
//	//for (int i = 0; i < _dims; i++)
//	//{
//	//	sizeImage[i] = imsizeX;
//	//	sizeImage[i] = imsizeX;
//	//}
//	//sizeImage[0] = dim[0];
//	//sizeImage[1] = dim[1];
//	cout << "image1 " << image1.width() << " " << image1.height() << endl;
//	cout << "image2 " << image2.width() << " " << image2.height() << endl;
//	cout << "sizeImage " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
//	//cout << "Input " << _inputV(30,30,0,0) << " " << _inputV(80, 80, 0, 0) << endl;
//	cout << " in the L1TVOpticalFlowNonlinear " <<endl;
//
//	/*float *inputV = 0;
//	float *inputY = 0;*/
//	//image  inputV = _inputV;
//	//image  inputY = _inputY;
//
//	int typeNorm = _norm;
//	//if (nrhs > 5)
//	//{
//	//	typeNorm = (int)mxGetScalar(prhs[5]);
//	//}
//
//
//	//if (nrhs > 6)
//	//{
//	//	inputV = mxGetPr(prhs[6]);
//	//	//mexPrintf("Given input vector field\n");
//	//}
//
//	//if (nrhs > 7)
//	//{
//	//	inputY = mxGetPr(prhs[7]);
//	//}
//
//	//inputV = _inputV;
//	//inputY = _inputY;
//
//	float stepsize[3] = { 1.0f, 1.0f, 1.0f };
//	//float stepsize[3] = { _stepsize[0], _stepsize[1], _stepsize[2] };
//	//if (nrhs > 8)
//	//{
//	//	double *tmpStepsize = mxGetPr(prhs[8]);
//
//	//	stepsize[0] = (float)tmpStepsize[0];
//	//	stepsize[1] = (float)tmpStepsize[1];
//	//	stepsize[2] = (float)tmpStepsize[2];
//	//}
//	float stepsizeD[3] = { 1.0f / stepsize[0], 1.0f / stepsize[1], 1.0f / stepsize[2] };
//
//	int numberOfWarps = _numberOfWarps;
//	//if (nrhs > 10)
//	//{
//	//	numberOfWarps = (int)mxGetScalar(prhs[10]);
//	//}
//
//	float huberEpsilon = 0.01f;
//	//if (nrhs > 11)
//	//{
//	//	huberEpsilon = (float)mxGetScalar(prhs[11]);
//	//}
//
//	int gradientConstancy = 0;
//	//if (nrhs > 12)
//	//{
//	//	gradientConstancy = (int)mxGetScalar(prhs[12]);
//	//}
//	//
//	const int nPx = (int)(dim[0] * dim[1]);
//
//	//const size_t sizeY[2] = { 7 * nPx, 1 };
//	const size_t sizeY[2] = { 5 * nPx, 1 };
//
//	// Output v1
//	//plhs[0] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
//	//double *Outv1 = mxGetPr(plhs[0]);
//
//	//// Output v2
//	//plhs[1] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
//	//double *Outv2 = mxGetPr(plhs[1]);
//
//	//// Output  Y
//	//plhs[2] = mxCreateNumericArray(2, sizeY, mxDOUBLE_CLASS, mxREAL);
//	//double *YOut = mxGetPr(plhs[2]);
//
//	float *Outv1 = new float[nPx];
//	float *Outv2 = new float[nPx];
//	//float *YOut = new float[nPx];
//
//	float* v1 = new float[nPx];
//	float* v2 = new float[nPx];
//
//	float* v1Old = new float[nPx];
//	float* v2Old = new float[nPx];
//
//	float* image1f = new float[nPx];
//	float* image2f = new float[nPx];
//
//	float* ux = new float[nPx];
//	float* uy = new float[nPx];
//	float* ut = new float[nPx];
////	image uXYZ(dim[0], dim[1], dim[2], 2, 0);
//	/*float* uxx = new float[nPx];
//	float* uxy = new float[nPx];
//	float* uyx = new float[nPx];
//	float* uyy = new float[nPx];*/
//
////	float* uxt = new float[nPx];
//	//float* uyt = new float[nPx];
//
//	float* y11 = new float[nPx];
//	float* y12 = new float[nPx];
//	float* y21 = new float[nPx];
//	float* y22 = new float[nPx];
//	float* y3 = new float[nPx];
//	/*image Y1(dim[0], dim[1], dim[2], 2, 0);
//	image Y2(dim[0], dim[1], dim[2], 2, 0);
//	image Y3(dim[0], dim[1], dim[2], 1, 0);*/
//	//float* y4 = new float[nPx];
//	//float* y5 = new float[nPx];
//
//	float* y11Old = new float[nPx];
//	float* y12Old = new float[nPx];
//	float* y21Old = new float[nPx];
//	float* y22Old = new float[nPx];
//	float* y3Old = new float[nPx];
//	//image Y1old(dim[0], dim[1], dim[2], 2, 0);
//	//image Y2old(dim[0], dim[1], dim[2], 2, 0);
//	//image Y3old(dim[0], dim[1], dim[2], 1, 0);
//	//float* y4Old = new float[nPx];
//	//float* y5Old = new float[nPx];
//
//	float* Kty1 = new float[nPx];
//	float* Kty2 = new float[nPx];
//
//	float* Kty1Old = new float[nPx];
//	float* Kty2Old = new float[nPx];
//
//	float* Kx11 = new float[nPx];
//	float* Kx12 = new float[nPx];
//	float* Kx21 = new float[nPx];
//	float* Kx22 = new float[nPx];
//	float* Kx3 = new float[nPx];
//	/*image KX1(dim[0], dim[1], dim[2], 2, 0);
//	image KX2(dim[0], dim[1], dim[2], 2, 0);
//	image KX3(dim[0], dim[1], dim[2], 1, 0);*/
//	//float* Kx4 = new float[nPx];
//	//float* Kx5 = new float[nPx];
//
//	float* Kx11Old = new float[nPx];
//	float* Kx12Old = new float[nPx];
//	float* Kx21Old = new float[nPx];
//	float* Kx22Old = new float[nPx];
//	float* Kx3Old = new float[nPx];
//
//	//image KX1old(dim[0], dim[1], dim[2], 2, 0);
//	//image KX2old(dim[0], dim[1], dim[2], 2, 0);
//	//image KX3old(dim[0], dim[1], dim[2], 1, 0);
//	//float* Kx4Old = new float[nPx];
//	//float* Kx5Old = new float[nPx];
//
//	float sigma1 = myMin(stepsize[1] / 2.0f, stepsize[2] / 2.0f);
//	float* sigma2 = new float[nPx];
//
//	//for gradient constancy
//	//float* sigma3 = new float[nPx];
//	//float* sigma4 = new float[nPx];
//
//	float* tau1 = new float[nPx];
//	float* tau2 = new float[nPx];
//
//	int * tableI = new int[nPx];
//	int * tableJ = new int[nPx];
//
//	//Huber Factor
//
//	const float huberFactor = 1.0f / (1.0f + sigma1* huberEpsilon / lambda);
//
//	//residuals
//	float p = 0.0f;
//	float d = 0.0f;
//
//	#pragma omp parallel for
//	for (int j = 0; j < dim[1]; ++j)
//	{
//		for (int i = 0; i < dim[0]; ++i)
//		{
//			int tmpIndex = index2DtoLinear(dim, i, j);
//
//			tableI[tmpIndex] = i;
//			tableJ[tmpIndex] = j;
//
//			//image1f[tmpIndex] = (float)u1[tmpIndex];
//			//image2f[tmpIndex] = (float)u2[tmpIndex];
//
//			/*	if (nrhs > 6)
//			{
//			v1[tmpIndex] = (float)inputV[tmpIndex];
//			v2[tmpIndex] = (float)inputV[nPx + tmpIndex];
//			}
//			else
//			{
//			v1[tmpIndex] = 0.0f;
//			v2[tmpIndex] = 0.0f;
//			}
//			*/
//			/*v1[tmpIndex] = (float)inputV[tmpIndex];
//			v2[tmpIndex] = (float)inputV[nPx + tmpIndex];*/
//
//	/*		v1[tmpIndex] = (float)inputV(tmpIndex);
//			v2[tmpIndex] = (float)inputV(nPx + tmpIndex);*/
//			v1[tmpIndex] = (float)_inputV(i, j, 0, 0);
//			v2[tmpIndex] = (float)_inputV(i, j, 0, 1);
//			//cout << "v1[tmpIndex]:	" << v1[tmpIndex] << " " << v2[tmpIndex] << endl;
//			//_inputY1(i, j, 0, 0);
//			//v1[tmpIndex] = 0.0f;
//			//v2[tmpIndex] = 0.0f;
//
//			Kty1[tmpIndex] = 0.0f;
//			Kty2[tmpIndex] = 0.0f;
//
//			//if (nrhs > 7)
//			//{
//			//	y11[tmpIndex] = (float)inputY[0 * nPx + tmpIndex];
//			//	y12[tmpIndex] = (float)inputY[1 * nPx + tmpIndex];
//			//	y21[tmpIndex] = (float)inputY[2 * nPx + tmpIndex];
//			//	y22[tmpIndex] = (float)inputY[3 * nPx + tmpIndex];
//			//             y3[tmpIndex] = (float)inputY[4 * nPx + tmpIndex];
//			//             
//			//             //for gradient constancy
//			//             y4[tmpIndex] = (float)inputY[5 * nPx + tmpIndex];
//			//             y5[tmpIndex] = (float)inputY[6 * nPx + tmpIndex];
//			//}
//			//else
//			//{
//			//	y11[tmpIndex] = 0.0f;
//			//	y12[tmpIndex] = 0.0f;
//			//	y21[tmpIndex] = 0.0f;
//			//	y22[tmpIndex] = 0.0f;
//			//             y3[tmpIndex] = 0.0f;
//			//             
//			//             //for gradient constancy
//			//             y4[tmpIndex] = 0.0f;
//			//             y5[tmpIndex] = 0.0f;
//			//}
//
//		/*	y11[tmpIndex] = (float)_inputY(0 * nPx + tmpIndex);
//			y12[tmpIndex] = (float)_inputY(1 * nPx + tmpIndex);
//			y21[tmpIndex] = (float)_inputY(2 * nPx + tmpIndex);
//			y22[tmpIndex] = (float)_inputY(3 * nPx + tmpIndex);
//			y3[tmpIndex] = (float)inputY(4 * nPx + tmpIndex);
//			*/
//			y11[tmpIndex] = (float)_inputY1(i, j, 0, 0);
//			y12[tmpIndex] = (float)_inputY1(i, j, 0, 1);
//			y21[tmpIndex] = (float)_inputY2(i, j, 0, 0);
//			y22[tmpIndex] = (float)_inputY2(i, j, 0, 1);
//			y3[tmpIndex] = (float)_inputY3(i, j, 0, 0);
//
//			//for gradient constancy
//			//	y4[tmpIndex] = (float)inputY(5 * nPx + tmpIndex);
//			//	y5[tmpIndex] = (float)inputY(6 * nPx + tmpIndex);
//
//
//
//			Kx11[tmpIndex] = 0.0f;
//			Kx12[tmpIndex] = 0.0f;
//			Kx21[tmpIndex] = 0.0f;
//			Kx22[tmpIndex] = 0.0f;
//			Kx3[tmpIndex] = 0.0f;
//
//			////for gradient constancy
//			//Kx4[tmpIndex] = 0.0f;
//			//Kx5[tmpIndex] = 0.0f;
//		}
//	}
//	clock_t begin = clock();
//	//do k warpings
//	for (int k = 0; k < numberOfWarps; ++k)
//	{
//	
//
//		// output??  u2, u21x,u21y, u22x,u22y
//		//		doWarp(  ut, ux, uy,uxt,uyt,uxx,uxy,uyx,uyy, sizeImage);
//		// output: ux, uy,  uxx,uxy,uyx,uyy
//		//  uxt,uyt : maybe for gradient consistency,
//		//ut,  
//
//		//doWarp(image1f, image2f, v1, v2, ut, ux, uy,uxt,uyt,uxx,uxy,uyx,uyy, sizeImage);
//		//doWarp(image1f, image2f, v1, v2, ut, ux, uy, sizeImage);
//
//
//		//////image test1(v1, sizeImage[0], sizeImage[1], 1, 1);
//		//////image test2(v2, sizeImage[0], sizeImage[1], 1, 1);
//
//		real	m_aabb[] = {-584,584,- 388,388,-10,10};
//		real dt = 1;
//		image _warpI2 = advection::advect(dim, m_aabb, _inputV, image2, -dt);
//			imagelist warpI2_xy = _warpI2.get_gradient("xy",0);
//		//	image im_ux = warpI2_xy[0];
//		//	image im_uy = warpI2_xy[1];
//			#pragma omp parallel for
//			for (int j = 0; j < dim[1]; ++j)
//			{
//				for (int i = 0; i < dim[0]; ++i)
//				{
//					int tmpIndex = index2DtoLinear(dim, i, j);
//					//ux[tmpIndex] = im_ux(tmpIndex);
//					//uy[tmpIndex] = im_uy(tmpIndex);
//					////cout << ux[tmpIndex] << " " << uy[tmpIndex] << endl;
//					//ut[tmpIndex] = warpIm2(tmpIndex) - image1(tmpIndex) - ux[tmpIndex]*v1[tmpIndex] - uy[tmpIndex]*v2[tmpIndex];
//					//uXYZ(i,j,0,0)= warpI2_xy[0](i, j);
//					//uXYZ(i, j, 0, 1) = warpI2_xy[1](i, j);
//					ux[tmpIndex] = warpI2_xy[0](i, j);
//					uy[tmpIndex] = warpI2_xy[1](i, j);
//					//cout << ux[tmpIndex] << " " << uy[tmpIndex] << endl;
//					//ut[tmpIndex] = _warpI2(i, j) - image1(i,j) - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex];
//					ut[tmpIndex] = _warpI2(i, j) - image1(i, j) - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex];
//					//ut[tmpIndex] = _warpI2(i, j) - image1(i, j) - uXYZ(i, j, 0, 0) * v1[tmpIndex] - uXYZ(i, j, 0, 1) * v2[tmpIndex];
//
//				}
//			}
//
//			cout << " in the L1TVOpticalFlowNonlinear  2 " << endl;
//
//
//
//
//
//		#pragma omp parallel for
//		for (int i = 0; i < nPx; ++i)
//		{
//			//adjust step sizes
//			
//			tau1[i] = 4.0f / std::min(stepsize[1], stepsize[2]) + myAbs(ux[i]);
//			tau2[i] = 4.0f / std::min(stepsize[1], stepsize[2]) + myAbs(uy[i]);
//			//tau1[i] = 4.0f / std::min(stepsize[1], stepsize[2]) + myAbs(uXYZ(i));
//			//tau2[i] = 4.0f / std::min(stepsize[1], stepsize[2]) + myAbs(uXYZ(i + nPx));
//			/*    if (gradientConstancy>0)
//			{
//			tau1[i] += std::abs(uxx[i]) + std::abs(uyx[i]);
//			tau2[i] += std::abs(uxy[i]) + std::abs(uyy[i]);
//			}*/
//			//sigma2[i] = std::abs(ux[i]) + std::abs(uXYZ(i + nPx));
//			sigma2[i] = std::abs(ux[i]) + std::abs(uy[i]);
//			//sigma2[i] = std::abs(uXYZ(i)) + std::abs(uXYZ(i+ nPx));
//
//			//for gradient constancy
//			/*         sigma3[i] = std::abs(uxx[i]) + std::abs(uxy[i]);
//			sigma4[i] = std::abs(uyx[i]) + std::abs(uyy[i]);
//			*/
//			tau1[i] = 1.0f / tau1[i];
//			tau2[i] = 1.0f / tau2[i];
//			sigma2[i] = 1.0f / sigma2[i];
//
//			//for gradient constancy
//			//sigma3[i] = 1.0f / sigma3[i];
//			//sigma4[i] = 1.0f / sigma4[i];
//		}
//
//
//		int iterations = 0;
//		float err = 1.0f;
//
//		while (err > tol && iterations <= maxIterations)
//		{
//			++iterations;
//
//			if (iterations % 50 == 0)
//			{
//				p = 0.0f;
//				d = 0.0f;
//			}
//
//			//cout << " in the L1TVOpticalFlowNonlinear  3  //primal step" << endl;
//
//			//primal step
//			#pragma omp parallel for
//			for (int j = 0; j < dim[1]; ++j)
//			{
//				for (int i = 0; i < dim[0]; ++i)
//				{
//					int tmpIndex = index2DtoLinear(dim, i, j);
//
//					/*imagelist Y11grad = _inputY1.get_shared_channel(0).get_gradient("xy", -1);
//					imagelist Y12grad = _inputY1.get_shared_channel(1).get_gradient("xy", -1);
//					imagelist Y21grad = _inputY2.get_shared_channel(0).get_gradient("xy", -1);
//					imagelist Y22grad = _inputY2.get_shared_channel(1).get_gradient("xy", -1);*/
//					// we have some issues here.
//					/*#pragma omp parallel for
//					for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
//					{*/
//						if (iterations % 50 == 0)
//						{
//							v1Old[tmpIndex] = v1[tmpIndex];
//							v2Old[tmpIndex] = v2[tmpIndex];
//
//							Kty1Old[tmpIndex] = Kty1[tmpIndex];
//							Kty2Old[tmpIndex] = Kty2[tmpIndex];
//						}
//
//
//						//transpose equals -div  // m is mb  not choose 1
//					/*	Kty1[tmpIndex] = -stepsizeD[1] * Y11grad[0](tmpIndex)
//										 - stepsizeD[2] * Y12grad[1](tmpIndex) + uXYZ(tmpIndex) * _inputY3(tmpIndex);
//
//						Kty2[tmpIndex] = -stepsizeD[1] * Y21grad[0](tmpIndex)
//										 - stepsizeD[2] * Y22grad[1](tmpIndex) + uXYZ(tmpIndex + nPx)  * _inputY3(tmpIndex);*/
//
//										 //transpose equals -div  // m is mb  not choose 1
//						Kty1[tmpIndex] = -stepsizeD[1] * dxm(y11, dim, tableI[tmpIndex], tableJ[tmpIndex])
//							- stepsizeD[2] * dym(y12, dim, tableI[tmpIndex], tableJ[tmpIndex]) + ux[tmpIndex] * y3[tmpIndex];
//
//						Kty2[tmpIndex] = -stepsizeD[1] * dxm(y21, dim, tableI[tmpIndex], tableJ[tmpIndex])
//							- stepsizeD[2] * dym(y22, dim, tableI[tmpIndex], tableJ[tmpIndex]) + uy[tmpIndex] * y3[tmpIndex];
//
//						//if (gradientConstancy>0)  // add a prior... means 
//						//{
//						//    Kty1[i] += uxx[i] * y4[i] + uyx[i] * y5[i];
//						//    Kty2[i] += uxy[i] * y4[i] + uyy[i] * y5[i];
//						//}
//
//						v1[tmpIndex] -= tau1[tmpIndex] * Kty1[tmpIndex];
//						v2[tmpIndex] -= tau2[tmpIndex] * Kty2[tmpIndex];
//						//_inputV(i,j,0,0) = v1[tmpIndex];// -= tau1[i] * Kty1[i];
//						//_inputV(i,j,0,1) = v2[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];
//						_inputV(tmpIndex) = v1[tmpIndex];// -= tau1[i] * Kty1[i];
//						_inputV(tmpIndex + nPx) = v2[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];
//
//						//_inputV(tmpIndex) -= tau1[tmpIndex] * Kty1[tmpIndex];
//						//_inputV(tmpIndex+ nPx) -= tau2[tmpIndex] * Kty2[tmpIndex];
//
//
//						if (iterations % 50 == 0)
//						{
//							//residuals
//							p += std::abs((v1Old[tmpIndex] - v1[tmpIndex]) / tau1[tmpIndex] - Kty1Old[tmpIndex] + Kty1[tmpIndex])
//								+ std::abs((v2Old[tmpIndex] - v2[tmpIndex]) / tau2[tmpIndex] - Kty2Old[tmpIndex] + Kty2[tmpIndex]);
//
//						}
//					}
//				}
//
//			//}		//cout << " in the L1TVOpticalFlowNonlinear  4  //dual step" << endl;
//
//			//dual step
//			//imagelist V1grad=_inputV.get_shared_channel(0).get_gradient("xy", 1);
//			//imagelist V2grad = _inputV.get_shared_channel(1).get_gradient("xy", 1);
//			/*imagelist V1grad = _inputV.get_shared_channel(0).get_gradient("xy", 1);
//			imagelist V2grad = _inputV.get_shared_channel(1).get_gradient("xy", 1);*/
//
//			//imagelist Vgrad = _inputV.get_gradient("xy", 1);
//			//image im1 = _inputV
//			#pragma omp parallel for reduction(+:d)
//			for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
//			{
//				float y11Tilde, y12Tilde, y21Tilde, y22Tilde;
//
//				if (iterations % 50 == 0)
//				{
//				
//					
//					/*Y1old(tmpIndex) = _inputY1(tmpIndex);
//					Y1old(tmpIndex+ nPx) = _inputY1(tmpIndex+ nPx);
//
//					Y2old(tmpIndex) = _inputY2(tmpIndex);
//					Y2old(tmpIndex + nPx) = _inputY2(tmpIndex + nPx);
//					Y3old(tmpIndex) = _inputY3(tmpIndex);*/
//					y11Old[tmpIndex] = y11[tmpIndex];
//					y12Old[tmpIndex] = y12[tmpIndex];
//					y21Old[tmpIndex] = y21[tmpIndex];
//					y22Old[tmpIndex] = y22[tmpIndex];
//
//					y3Old[tmpIndex] = y3[tmpIndex];
//				}
//
//				/*KX1old(tmpIndex) = KX1(tmpIndex);
//				KX1old(tmpIndex+ nPx) = KX1(tmpIndex+ nPx);
//				KX2old(tmpIndex) = KX2(tmpIndex);
//				KX2old(tmpIndex + nPx) = KX2(tmpIndex + nPx);
//				KX3old(tmpIndex) = KX3(tmpIndex);*/
//				Kx11Old[tmpIndex] = Kx11[tmpIndex];
//				Kx12Old[tmpIndex] = Kx12[tmpIndex];
//				Kx21Old[tmpIndex] = Kx21[tmpIndex];
//				Kx22Old[tmpIndex] = Kx22[tmpIndex];
//
//				Kx3Old[tmpIndex] = Kx3[tmpIndex];
//			/*	
//				KX1(tmpIndex) = stepsizeD[1] * V1grad[0](tmpIndex);
//				KX1(tmpIndex + nPx) = stepsizeD[2] * V1grad[1](tmpIndex);
//				KX2(tmpIndex) = stepsizeD[1] * V2grad[0](tmpIndex);
//				KX2(tmpIndex + nPx) = stepsizeD[2] * V2grad[1](tmpIndex);*/
//				/*Kx11[tmpIndex] = stepsizeD[1] * V1grad[0](tmpIndex);
//				Kx12[tmpIndex] = stepsizeD[2] * V1grad[1](tmpIndex);
//				Kx21[tmpIndex] = stepsizeD[1] * V2grad[0](tmpIndex);
//				Kx22[tmpIndex] = stepsizeD[2] * V2grad[1](tmpIndex);*/
//				Kx11[tmpIndex] = stepsizeD[1] * dxp(v1, dim, tableI[tmpIndex], tableJ[tmpIndex]);
//				Kx12[tmpIndex] = stepsizeD[2] * dyp(v1, dim, tableI[tmpIndex], tableJ[tmpIndex]);
//				Kx21[tmpIndex] = stepsizeD[1] * dxp(v2, dim, tableI[tmpIndex], tableJ[tmpIndex]);
//				Kx22[tmpIndex] = stepsizeD[2] * dyp(v2, dim, tableI[tmpIndex], tableJ[tmpIndex]);
//
//				Kx3[tmpIndex] = ux[tmpIndex] * v1[tmpIndex] + uy[tmpIndex] * v2[tmpIndex];  // warpped_u2=(ux, uy)
//
//				//KX3(tmpIndex) = uXYZ(tmpIndex) * v1[tmpIndex] + uXYZ(tmpIndex + nPx) * v2[tmpIndex];  // warpped_u2=(ux, uy)
//				//Kx3[tmpIndex] = uXYZ(tmpIndex) * v1[tmpIndex] + uXYZ(tmpIndex+ nPx) * v2[tmpIndex];  // warpped_u2=(ux, uy)
//
//				/*Kx11[tmpIndex] = stepsizeD[1] * dxp(_inputV.get_channel(0), dim, tableI[tmpIndex], tableJ[tmpIndex]);
//				Kx12[tmpIndex] = stepsizeD[2] * dyp(_inputV.get_channel(0), dim, tableI[tmpIndex], tableJ[tmpIndex]);
//				Kx21[tmpIndex] = stepsizeD[1] * dxp(_inputV.get_channel(1), dim, tableI[tmpIndex], tableJ[tmpIndex]);
//				Kx22[tmpIndex] = stepsizeD[2] * dyp(_inputV.get_channel(1), dim, tableI[tmpIndex], tableJ[tmpIndex]);*/
//
//				//Kx3[tmpIndex] = ux[tmpIndex] * _inputV(tmpIndex) + uy[tmpIndex] * _inputV(tmpIndex+ nPx);  // warpped_u2=(ux, uy)
//
//																							/*if (gradientConstancy>0)
//																							{
//																							Kx4Old[tmpIndex] = Kx4[tmpIndex];
//																							Kx5Old[tmpIndex] = Kx5[tmpIndex];
//
//																							Kx4[tmpIndex] = uxx[tmpIndex] * v1[tmpIndex] + uxy[tmpIndex] * v2[tmpIndex];
//																							Kx5[tmpIndex] = uyx[tmpIndex] * v1[tmpIndex] + uyy[tmpIndex] * v2[tmpIndex];
//
//																							y4Old[tmpIndex] = y4[tmpIndex];
//																							y5Old[tmpIndex] = y5[tmpIndex];
//
//																							y4[tmpIndex] = std::max(-1.0f,std::min(1.0f, y4[tmpIndex] + sigma3[tmpIndex]*(Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old[tmpIndex] + uxt[tmpIndex])));
//																							y5[tmpIndex] = std::max(-1.0f,std::min(1.0f, y5[tmpIndex] + sigma4[tmpIndex]*(Kx5[tmpIndex] + Kx5[tmpIndex] - Kx5Old[tmpIndex] + uyt[tmpIndex])));
//																							}*/
//
//				if (typeNorm == 4) // Huber
//				{
//					/*y11Tilde = (_inputY1(tmpIndex) + sigma1*(KX1(tmpIndex) + KX1(tmpIndex) - KX1old(tmpIndex))) * huberFactor;
//					y12Tilde = (_inputY1(tmpIndex + nPx) + sigma1*(KX1(tmpIndex + nPx) + KX1(tmpIndex + nPx) - KX1old(tmpIndex + nPx))) * huberFactor;
//					y21Tilde = (_inputY2(tmpIndex)+ sigma1*(KX2(tmpIndex) + KX2(tmpIndex) - KX2old(tmpIndex))) * huberFactor;
//					y22Tilde = (_inputY2(tmpIndex + nPx) + sigma1*(KX2(tmpIndex + nPx) + KX2(tmpIndex + nPx) - KX2old(tmpIndex + nPx))) * huberFactor;*/
//					y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex])) * huberFactor;
//					y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex])) * huberFactor;
//					y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex])) * huberFactor;
//					y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex])) * huberFactor;
//				}
//				else
//				{
//					/*y11Tilde = _inputY1(tmpIndex) + sigma1*(KX1(tmpIndex) + KX1(tmpIndex) - KX1old(tmpIndex));
//					y12Tilde = _inputY1(tmpIndex + nPx) + sigma1*(KX1(tmpIndex + nPx) + KX1(tmpIndex + nPx) - KX1old(tmpIndex + nPx));
//					y21Tilde = _inputY2(tmpIndex) + sigma1*(KX2(tmpIndex) + KX2(tmpIndex) - KX2old(tmpIndex));
//					y22Tilde = _inputY2(tmpIndex + nPx) + sigma1*(KX2(tmpIndex + nPx) + KX2(tmpIndex + nPx) - KX2old(tmpIndex + nPx)) ;*/
//
//					y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex]));
//					y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex]));
//					y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex]));
//					y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex]));
//				}
//
//				float divisor1 = std::max(1.0f, std::sqrt(y11Tilde*y11Tilde + y12Tilde*y12Tilde) / lambda);
//				float divisor2 = std::max(1.0f, std::sqrt(y21Tilde*y21Tilde + y22Tilde*y22Tilde) / lambda);
//
//				/*_inputY1(tmpIndex) = y11Tilde / divisor1;
//				_inputY1(tmpIndex+ nPx) = y12Tilde / divisor1;
//				_inputY2(tmpIndex) = y21Tilde / divisor2;
//				_inputY2(tmpIndex+ nPx) = y22Tilde / divisor2;*/
//				y11[tmpIndex] = y11Tilde / divisor1;
//				y12[tmpIndex] = y12Tilde / divisor1;
//				y21[tmpIndex] = y21Tilde / divisor2;
//				y22[tmpIndex] = y22Tilde / divisor2;
//
//
//
//				//_inputY3(tmpIndex) = std::max(-1.0f, std::min(1.0f, _inputY3(tmpIndex) + sigma2[tmpIndex] * (KX3(tmpIndex) + KX3(tmpIndex) - KX3old(tmpIndex) + ut[tmpIndex])));
//
//				// where to get the u_t, 
//				y3[tmpIndex] = std::max(-1.0f, std::min(1.0f, y3[tmpIndex] + sigma2[tmpIndex] * (Kx3[tmpIndex] + Kx3[tmpIndex] - Kx3Old[tmpIndex] + ut[tmpIndex])));
//			
//				if (iterations % 50 == 0)
//				{
//
//					//d += std::abs((Y1old(tmpIndex) - _inputY1(tmpIndex)) / sigma1 - KX1old(tmpIndex) + KX1(tmpIndex)) +
//					//	std::abs((Y1old(tmpIndex + nPx) - _inputY1(tmpIndex + nPx)) / sigma1 - KX1old(tmpIndex + nPx) + KX1(tmpIndex + nPx)) +
//					//	std::abs((Y2old(tmpIndex) - _inputY2(tmpIndex)) / sigma1 - KX2old(tmpIndex) + KX2(tmpIndex)) +
//					//	std::abs((Y2old(tmpIndex + nPx) - _inputY2(tmpIndex + nPx)) / sigma1 - KX2old(tmpIndex + nPx) + KX2(tmpIndex + nPx)) +
//					//	std::abs((Y3old(tmpIndex) - _inputY3(tmpIndex)) / sigma2[tmpIndex] - KX3old(tmpIndex) + KX3(tmpIndex));
//
//					d += std::abs((y11Old[tmpIndex] - y11[tmpIndex]) / sigma1 - Kx11Old[tmpIndex] + Kx11[tmpIndex]) +
//						std::abs((y12Old[tmpIndex] - y12[tmpIndex]) / sigma1 - Kx12Old[tmpIndex] + Kx12[tmpIndex]) +
//						std::abs((y21Old[tmpIndex] - y21[tmpIndex]) / sigma1 - Kx21Old[tmpIndex] + Kx21[tmpIndex]) +
//						std::abs((y22Old[tmpIndex] - y22[tmpIndex]) / sigma1 - Kx22Old[tmpIndex] + Kx22[tmpIndex]) +
//						std::abs((y3Old[tmpIndex] - y3[tmpIndex]) / sigma2[tmpIndex] - Kx3Old[tmpIndex] + Kx3[tmpIndex]);
//
//					/*if (gradientConstancy>0)
//					{
//						d += std::abs((y4Old[tmpIndex] - y4[tmpIndex]) / sigma3[tmpIndex] - Kx4Old[tmpIndex] + Kx4[tmpIndex]) +
//							std::abs((y5Old[tmpIndex] - y5[tmpIndex]) / sigma4[tmpIndex] - Kx5Old[tmpIndex] + Kx5[tmpIndex]);
//					}*/
//				}
//			}
//
//			if (iterations % 50 == 0)
//			{
//				err = (d*d + p*p) / (float)nPx;
//			}
//
//			if (iterations % 1000 == 0)
//			{
//
//				cout << "Iteration: " << iterations << " Residual " << err << endl;
//
//				//mexPrintf("Iteration %d,Residual %e\n", iterations, err);
//				//mexEvalString("drawnow;");
//			}
//		}
//	}
//	//cout << " in the L1TVOpticalFlowNonlinear 5 " << endl;
//
//	//write output
//	//#pragma omp parallel for
//	//for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
//	//{
//
//	// end of warps 
//
//	#pragma omp parallel for
//	for (int j = 0; j < dim[1]; ++j)
//	{
//		for (int i = 0; i < dim[0]; ++i)
//		{
//			int tmpIndex = index2DtoLinear(dim, i, j);
//			/*
//			YOut[tmpIndex + 0 * nPx] = (double)y11[tmpIndex];
//			YOut[tmpIndex + 1 * nPx] = (double)y12[tmpIndex];
//			YOut[tmpIndex + 2 * nPx] = (double)y21[tmpIndex];
//			YOut[tmpIndex + 3 * nPx] = (double)y22[tmpIndex];
//			YOut[tmpIndex + 4 * nPx] = (double)y3[tmpIndex];
//
//			YOut[tmpIndex + 5 * nPx] = (double)y4[tmpIndex];
//			YOut[tmpIndex + 6 * nPx] = (double)y5[tmpIndex];
//
//			Outv1[tmpIndex] = (double)v1[tmpIndex];
//			Outv2[tmpIndex] = (double)v2[tmpIndex];*/
//
//			_inputY1(i, j, 0, 0) = (double)y11[tmpIndex];
//			_inputY1(i, j, 0, 1) = (double)y12[tmpIndex];
//			_inputY2(i, j, 0, 0) = (double)y21[tmpIndex];
//			_inputY2(i, j, 0, 1) = (double)y22[tmpIndex];
//			_inputY3(i, j, 0, 0) = (double)y3[tmpIndex];
//
//			//	_inputY(tmpIndex + 5 * nPx) = (double)y4[tmpIndex];
//			//	_inputY(tmpIndex + 6 * nPx) = (double)y5[tmpIndex];
//
//			//	_inputV[tmpIndex] += (double)v1[tmpIndex];
//			//	_inputV[tmpIndex +  nPx] += (double)v2[tmpIndex];
//			_inputV(i, j, 0, 0) = (float)v1[tmpIndex];
//			_inputV(i, j, 0, 1) = (float)v2[tmpIndex];
//			//	cout << _inputV(tmpIndex) << " " << _inputV(tmpIndex + nPx) << endl;
//
//			Outv1[tmpIndex] = _inputV(i, j, 0, 0);
//			Outv2[tmpIndex] = _inputV(i, j, 0, 1);
//			//	_inputV(tmpIndex)
//			//cout << Outv1[tmpIndex] << " ";
//		}
//
//	}
//	cout << " in the L1TVOpticalFlowNonlinear  5 " << endl;
//
//	clock_t end = clock();
//	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	cout << "=========Total time is :" << elapsed_secs<< " s ======"<<endl;
//	int num1 = _inputV.width();
//	int num2 = _inputV.height();
//	//itoa(num1, snum1, 10);
//	//itoa(num2, snum2, 10);
//	//str = snum1+"h" + snum2+"-1";
//	//str+="h"+
//	char Mybuf[10];
//	sprintf(Mybuf, "%d", num1);
//	string Myb = Mybuf;
//
//	char Mybuf2[10];
//	sprintf(Mybuf2, "%d", num2);
//	string Myb2 = Mybuf2;
//
//	string _output;// = pOutput;
//	_output += "/";
//	string Mywritefile1 = Myb + "-" + Myb2 + "-1.txt";
//	string Mywritefile2 = Myb + "-" + Myb2 + "-2.txt";
//
//	char* MyOutputfile1 = _strdup(Mywritefile1.c_str());
//	char* MyOutputfile2 = _strdup(Mywritefile2.c_str());
//	//strdup(
//	FILE *fp;
//	//fp = fopen("test1.txt", "wb");
//	fp = fopen(MyOutputfile1, "wb");
//
//	/*if (fp == null)
//	return;*/
//	//char x[10] = "ABCDEFGHIJ";
//	fwrite(Outv1, sizeof(float), dim[0] * dim[1], fp);
//	fclose(fp);
//
//
//	FILE *fp2;
//	fp2 = fopen(MyOutputfile2, "wb");
//	/*if (fp == null)
//	return;*/
//	//char x[10] = "ABCDEFGHIJ";
//	fwrite(Outv2, sizeof(float), dim[0] * dim[1], fp2);
//	fclose(fp2);
//
//	cout << "Testing the application end part" << endl;
//
//
//	//Sleep(1000);
//
//
//	//image warpField2(sizeImage[0], sizeImage[1], 1, 2);
//
//	////cimg_forXY(warpField, x, y) {
//	//for (int j = 0; j < sizeImage[1]; j++)
//	//{
//	//	for (int i = 0; i < sizeImage[0]; i++)
//	//	{
//	//		int tmpIndex = index2DtoLinear(sizeImage, i, j);
//	//		//const float
//	//		//u = x + x0 - dw2, v = y + y0 - dh2;
//
//	//		warpField2(i, j, 0) = -Outv1[tmpIndex];
//	//		warpField2(i, j, 1) = -Outv2[tmpIndex];
//	//	}
//	//}
//	//image ta = image2;
//	//real _max = image1.max();
//	//real _min= image1.min();
//	//image warpIm22 = ta.get_warp(warpField2, 1, 2, 0).normalize(_min, _max);
//
//
//	//image test2(Outv1, sizeImage[0], sizeImage[1], 1, 1);
//	//image test3(Outv2, sizeImage[0], sizeImage[1], 1, 1);
//	////test2 = -test2;
//	////test3 = -test3;
//	////cout << i2(0) << " " << i2(4) << " " << i2(8) << endl;
//	//image testwarpfield2 = test2.get_append(test3,'c');
//	////testwarpfield2.insert(test2);
//	////testwarpfield2.insert(test3);
//	////image testfiled = testwarpfield2.get_append();
//	//image tt = image2;
//	//image warpIm2t = tt.get_warp(testwarpfield2, 1, 2,1);
//
//	////image testwarpfield = image(v1).get_append(testymap, 'c');
//	/*
//	imagelist testresultlist(image2, image1, warpIm22);
//
//	Sleep(500);
//
//	testresultlist.display();
//	Sleep(500);
//	*/
//
//	/*	CImgList<float> teste(abs(image2 - warpIm2), abs(image2 - warpIm2));
//	teste.display();*/
//
//
//
//	delete[] tableI;
//	delete[] tableJ;
//
//	delete[] image1f;
//	delete[] image2f;
//
//	delete[] sigma2;
//	//delete[] sigma3;
//	//delete[] sigma4;
//	delete[] tau1;
//	delete[] tau2;
//
//	delete[] v1;
//	delete[] v2;
//	delete[] v1Old;
//	delete[] v2Old;
//
//	//delete[] ux;
//	//delete[] uy;
//	delete[] ut;
//
//	//delete[] uxx;
//	//delete[] uxy;
//	//delete[] uyx;
//	//delete[] uyy;
//
//	//delete[] uxt;
//	//delete[] uyt;
//
//	delete[] y11;
//	delete[] y12;
//	delete[] y21;
//	delete[] y22;
//	delete[] y3;
//	//delete[] y4;
//	//delete[] y5;
//
//	delete[] y11Old;
//	delete[] y12Old;
//	delete[] y21Old;
//	delete[] y22Old;
//	delete[] y3Old;
//	//delete[] y4Old;
//	//delete[] y5Old;
//
//	delete[] Kty1;
//	delete[] Kty2;
//
//	delete[] Kty1Old;
//	delete[] Kty2Old;
//
//	//delete[] Kx11;
//	//delete[] Kx12;
//	//delete[] Kx21;
//	//delete[] Kx22;
//	//delete[] Kx3;
//	//delete[] Kx4;
//	//delete[] Kx5;
//
//	//delete[] Kx11Old;
//	//delete[] Kx12Old;
//	//delete[] Kx21Old;
//	//delete[] Kx22Old;
//	//delete[] Kx3Old;
//	//delete[] Kx4Old;
//	//delete[] Kx5Old;
//
//	//return ;
//}
//
//#elif DIMENSION_SIZE ==3 // for 3d case
//void L1TVOpticalFlowNonlinear_MultiScale3D(const int *dim, const image &image1, const image &image2, image &_inputV,
//	image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4, float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps,
//	float eta, float sigma, int scales)
//{
//	// the last finest level should use blurred one or not.
//	image I1_blur = image1.get_blur(sigma);
//	image I2_blur = image2.get_blur(sigma);
//
//	//image tempI1 = I1_blur;
//	//image tempI2 = I2_blur;
//	//image tempV = _inputV;
//	//image tempY1 = _inputY1;
//	//image tempY2 = _inputY2;
//	//image tempY3 = _inputY3;
//	//image tempY4 = _inputY4;
//
//	//image I3_blur = image3.get_blur(sigma);
//	//float **Im1s;
//	//float **Im2s;
//	int size = image1.width()*image1.height()*image1.depth();
//	//Im1s[0] = new float[size];
//	//Im2s[0] = new float[size];
//	imagelist Im1s;
//	imagelist Im2s;
//	imagelist inputVs;
//	imagelist inputY1s;
//	imagelist inputY2s;
//	imagelist inputY3s;
//	imagelist inputY4s;
//
//	imagelist dims;
//	int twidth = I1_blur.width();
//	int theight = I1_blur.height();
//	int tdepth = I1_blur.depth();
//	float t_eta = 1.0f;
//	for (int i = 0; i < scales; i++)
//	{
//		Im1s.insert(I1_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//		//Im1s.insert()
//		//image timg = tempI1.get_resize(tempI1.width()*eta, tempI1.height()*eta, tempI1.depth()*eta, 1, 5);
//		//tempI1 = timg;
//		Im2s.insert(I2_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//		//Im2s.insert(tempI2);
//		//image timg2 = tempI2.get_resize(tempI2.width()*eta, tempI2.height()*eta, tempI2.depth()*eta, 1, 5);
//		//tempI2 = timg2;
//		inputVs.insert(_inputV.get_resize(twidth*t_eta, theight*t_eta, tdepth, 3, 5));
//
//		/*inputVs.insert(tempV);
//		image timg3 = tempV.get_resize(tempV.width()*eta, tempV.height()*eta, tempV.depth()*eta, 3, 5);
//		tempV = timg3;*/
//		inputY1s.insert(_inputY1.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, 5));
//		/*inputY1s.insert(tempY1);
//		image timg4 = tempY1.get_resize(tempY1.width()*eta, tempY1.height()*eta, tempY1.depth()*eta, 3, 5);
//		tempY1 = timg4;*/
//
//		inputY2s.insert(_inputY2.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, 5));
//	/*	inputY2s.insert(tempY2);
//		image timg5 = tempY2.get_resize(tempY2.width()*eta, tempY2.height()*eta, tempY2.depth()*eta, 3, 5);
//		tempY2 = timg5;*/
//
//		inputY3s.insert(_inputY3.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, 5));
//
//		//inputY3s.insert(tempY3);
//		//image timg6 = tempY3.get_resize(tempY3.width()*eta, tempY3.height()*eta, tempY3.depth()*eta, 3, 5);
//		//tempY3 = timg6;
//
//		inputY4s.insert(_inputY4.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//		/*inputY4s.insert(tempY4);
//		image timg7 = tempY4.get_resize(tempY4.width()*eta, tempY4.height()*eta, tempY4.depth()*eta, 1, 5);
//		tempY3 = timg6;*/
//
//		t_eta = t_eta*eta;
//	}
//	Im1s.reverse();
//	Im2s.reverse();
//	inputVs.reverse();
//	//inputYs.reverse();
//	inputY1s.reverse();
//	inputY2s.reverse();
//	inputY3s.reverse();
//	inputY4s.reverse();
//
//	////for (CImgList<>::iterator it = list.begin(); it<list.end(); ++it) (*it).mirror('x');
//	for (int i = 0; i < inputVs.size(); i++)
//	{
//		cout << "At levels " << inputVs.size() - i << endl;
//		int _dim[] = { Im1s(i).width(),Im1s(i).height(), Im1s(i).depth()};
//		//L1TVOpticalFlowNonlinear(_dim, Im1s(i), Im2s(i), inputVs(i), inputYs(i), _tol, _lambda, _maxIterations, _norm,
//		//									_numberOfWarps);
//		L1TVOpticalFlowNonlinear3D(_dim, Im1s(i), Im2s(i), inputVs(i), inputY1s(i), inputY2s(i), inputY3s(i), inputY4s(i), _tol, _lambda, _maxIterations, _norm,
//			_numberOfWarps);
//		if (i < inputVs.size() - 1)
//		{
//			inputVs(i + 1) = inputVs(i).get_resize(inputVs(i + 1).width(), inputVs(i + 1).height(), inputVs(i + 1).depth(), 3, 5);
//			//inputVs(i + 1) = inputVs(i).get_resize(inputVs(i + 1).width(), inputVs(i + 1).height(), inputVs(i + 1).depth(), 2, 5);
//			//inputYs(i + 1) = inputYs(i).get_resize(inputYs(i + 1).width(), inputYs(i + 1).height(), inputYs(i + 1).depth(), 1, 5);
//			inputY1s(i + 1) = inputY1s(i).get_resize(inputY1s(i + 1).width(), inputY1s(i + 1).height(), inputY1s(i + 1).depth(), 3, 5);
//			inputY2s(i + 1) = inputY2s(i).get_resize(inputY2s(i + 1).width(), inputY2s(i + 1).height(), inputY2s(i + 1).depth(), 3, 5);
//			inputY3s(i + 1) = inputY3s(i).get_resize(inputY3s(i + 1).width(), inputY3s(i + 1).height(), inputY3s(i + 1).depth(), 3, 5);
//			inputY4s(i + 1) = inputY4s(i).get_resize(inputY4s(i + 1).width(), inputY4s(i + 1).height(), inputY4s(i + 1).depth(), 1, 5);
//
//		}
//
//
//		//cout << inputVs.size()<<" "<< inputVs(i).width() << " " << inputVs(i).height() << endl;
//
//	}
//	////cout << Im1s.size() << endl;
//	////system("pause");
//
//}
//#else
//
//
//void L1TVOpticalFlowNonlinear3D(int *dim, const image & image1, const image & image2, image & _inputV, image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4,
//	float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps, real *_aabb, float huber)
//	//void L1TVOpticalFlowNonlinear3D(int *dim, const image & image1, const image & image2, image & _inputV1, image & _inputV2, image & _inputV3, 
//	//	image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4,
//	//	float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps, real *_aabb)
//{
//
//	//const int maxIterations = (int)mxGetScalar(prhs[4]);
//
//	//const size_t *sizeImage = mxGetDimensions(prhs[0]);
//
//	//float *u1 = image1;
//	//float *u2 = image2;
//	//image u1 = image1;
//	//image u2 = image2;
//
//
//	cout << "AABB: " << _aabb[0] << " " << _aabb[1] << " " << _aabb[2] << " " << _aabb[3] << " " << _aabb[4] << " " << _aabb[5] << " " << endl;
//	const float tol = _tol;
//	const float lambda = _lambda;
//
//	//mexPrintf("tol: %f, lambda: %f; ", tol, lambda);
//
//	const int maxIterations = _maxIterations;
//
//	//const size_t *sizeImage = mxGetDimensions(prhs[0]);
//	//int *sizeImage = new int[2];
//	//for (int i = 0; i < _dims; i++)
//	//{
//	//	sizeImage[i] = imsizeX;
//	//	sizeImage[i] = imsizeX;
//	//}
//	//sizeImage[0] = dim[0];
//	//sizeImage[1] = dim[1];
//	cout << "image1 " << image1.width() << " " << image1.height() << " " << image1.depth() << endl;
//	cout << "image2 " << image2.width() << " " << image2.height() << " " << image2.depth() << endl;
//	cout << "sizeImage " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
//	//cout << "Input " << _inputV(30,30,0,0) << " " << _inputV(80, 80, 0, 0) << endl;
//	cout << " in the L1TVOpticalFlowNonlinear " << endl;
//
//	/*float *inputV = 0;
//	float *inputY = 0;*/
//	//image  inputV = _inputV;
//	//image  inputY = _inputY;
//
//	int typeNorm = _norm;
//
//
//	float stepsize[4] = { 1.0f, 1.0f, 1.0f ,1.0f };
//	//float stepsize[3] = { _stepsize[0], _stepsize[1], _stepsize[2] };
//
//	float stepsizeD[4] = { 1.0f / stepsize[0],
//		1.0f / stepsize[1],
//		1.0f / stepsize[2],
//		1.0f / stepsize[3] };
//
//	int numberOfWarps = _numberOfWarps;
//	//if (nrhs > 10)
//	//{
//	//	numberOfWarps = (int)mxGetScalar(prhs[10]);
//	//}
//
//	//float huberEpsilon = 0.01f;
//	float huberEpsilon = huber;
//
//	//float huberEpsilon = 0.1f;
//	//if (nrhs > 11)
//	//{
//	//	huberEpsilon = (float)mxGetScalar(prhs[11]);
//	//}
//
//	int gradientConstancy = 0;
//	//if (nrhs > 12)
//	//{
//	//	gradientConstancy = (int)mxGetScalar(prhs[12]);
//	//}
//	//
//	const int nPx = (int)(dim[0] * dim[1] * dim[2]);
//	cout << "nPX" << nPx << endl;
//	//const size_t sizeY[2] = { 7 * nPx, 1 };
//	//const size_t sizeY[2] = { 5 * nPx, 1 };
//
//	// Output v1
//	//plhs[0] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
//	//double *Outv1 = mxGetPr(plhs[0]);
//
//	//// Output v2
//	//plhs[1] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
//	//double *Outv2 = mxGetPr(plhs[1]);
//
//	//// Output  Y
//	//plhs[2] = mxCreateNumericArray(2, sizeY, mxDOUBLE_CLASS, mxREAL);
//	//double *YOut = mxGetPr(plhs[2]);
//
//	//float *Outv1 = new float[nPx];
//	//float *Outv2 = new float[nPx];
//	//float *Outv3 = new float[nPx];
//	////float *YOut = new float[nPx];
//
//	float* v1 = new float[nPx];
//	float* v2 = new float[nPx];
//	float* v3 = new float[nPx];
//
//	float* v1Old = new float[nPx];
//	float* v2Old = new float[nPx];
//	float* v3Old = new float[nPx];
//
//	float* image1f = new float[nPx];
//	float* image2f = new float[nPx];
//
//	float* ux = new float[nPx];
//	float* uy = new float[nPx];
//	float* uz = new float[nPx];
//	float* ut = new float[nPx];
//	//	image uXYZ(dim[0], dim[1], dim[2], 2, 0);
//	/*float* uxx = new float[nPx];
//	float* uxy = new float[nPx];
//	float* uyx = new float[nPx];
//	float* uyy = new float[nPx];*/
//
//	//	float* uxt = new float[nPx];
//	//float* uyt = new float[nPx];
//
//	float* y11 = new float[nPx];
//	float* y12 = new float[nPx];
//	float* y13 = new float[nPx];
//
//	float* y21 = new float[nPx];
//	float* y22 = new float[nPx];
//	float* y23 = new float[nPx];
//
//	float* y31 = new float[nPx];
//	float* y32 = new float[nPx];
//	float* y33 = new float[nPx];
//	float* y4 = new float[nPx];
//	/*image Y1(dim[0], dim[1], dim[2], 2, 0);
//	image Y2(dim[0], dim[1], dim[2], 2, 0);
//	image Y3(dim[0], dim[1], dim[2], 1, 0);*/
//	//float* y4 = new float[nPx];
//	//float* y5 = new float[nPx];
//
//	float* y11Old = new float[nPx];
//	float* y12Old = new float[nPx];
//	float* y13Old = new float[nPx];
//
//	float* y21Old = new float[nPx];
//	float* y22Old = new float[nPx];
//	float* y23Old = new float[nPx];
//
//	float* y31Old = new float[nPx];
//	float* y32Old = new float[nPx];
//	float* y33Old = new float[nPx];
//	float* y4Old = new float[nPx];
//	//image Y1old(dim[0], dim[1], dim[2], 2, 0);
//	//image Y2old(dim[0], dim[1], dim[2], 2, 0);
//	//image Y3old(dim[0], dim[1], dim[2], 1, 0);
//	//float* y4Old = new float[nPx];
//	//float* y5Old = new float[nPx];
//
//	float* Kty1 = new float[nPx];
//	float* Kty2 = new float[nPx];
//	float* Kty3 = new float[nPx];
//
//	float* Kty1Old = new float[nPx];
//	float* Kty2Old = new float[nPx];
//	float* Kty3Old = new float[nPx];
//
//	float* Kx11 = new float[nPx];
//	float* Kx12 = new float[nPx];
//	float* Kx13 = new float[nPx];
//
//	float* Kx21 = new float[nPx];
//	float* Kx22 = new float[nPx];
//	float* Kx23 = new float[nPx];
//
//	float* Kx31 = new float[nPx];
//	float* Kx32 = new float[nPx];
//	float* Kx33 = new float[nPx];
//	float* Kx4 = new float[nPx];
//	/*image KX1(dim[0], dim[1], dim[2], 2, 0);
//	image KX2(dim[0], dim[1], dim[2], 2, 0);
//	image KX3(dim[0], dim[1], dim[2], 1, 0);*/
//	//float* Kx4 = new float[nPx];
//	//float* Kx5 = new float[nPx];
//
//	float* Kx11Old = new float[nPx];
//	float* Kx12Old = new float[nPx];
//	float* Kx13Old = new float[nPx];
//
//	float* Kx21Old = new float[nPx];
//	float* Kx22Old = new float[nPx];
//	float* Kx23Old = new float[nPx];
//
//	float* Kx31Old = new float[nPx];
//	float* Kx32Old = new float[nPx];
//	float* Kx33Old = new float[nPx];
//	float* Kx4Old = new float[nPx];
//
//	//image KX1old(dim[0], dim[1], dim[2], 2, 0);
//	//image KX2old(dim[0], dim[1], dim[2], 2, 0);
//	//image KX3old(dim[0], dim[1], dim[2], 1, 0);
//	//float* Kx4Old = new float[nPx];
//	//float* Kx5Old = new float[nPx];
//
//	//float sigma1 = myMin(stepsize[1] / 2.0f, stepsize[2] / 2.0f);
//
//	float sigma1 = myMin(stepsize[3] / 3.0f, myMin(stepsize[1] / 3.0f, stepsize[2] / 3.0f));
//	float* sigma2 = new float[nPx];
//
//	//for gradient constancy
//	//float* sigma3 = new float[nPx];
//	//float* sigma4 = new float[nPx];
//
//	float* tau1 = new float[nPx];
//	float* tau2 = new float[nPx];
//	float* tau3 = new float[nPx];
//
//	int * tableI = new int[nPx];
//	int * tableJ = new int[nPx];
//	int * tableK = new int[nPx];
//
//	//Huber Factor
//
//	const float huberFactor = 1.0f / (1.0f + sigma1* huberEpsilon / lambda);
//
//	//residuals
//	float p = 0.0f;
//	float d = 0.0f;
//
//	#pragma omp parallel for
//	for (int k = 0; k < dim[2]; ++k)
//	{
//		for (int j = 0; j < dim[1]; ++j)
//		{
//			for (int i = 0; i < dim[0]; ++i)
//			{
//				//int tmpIndex = index2DtoLinear(dim, i, j);
//				int tmpIndex = index3DtoLinear(dim, i, j, k);
//
//				tableI[tmpIndex] = i;
//				tableJ[tmpIndex] = j;
//				tableK[tmpIndex] = k;
//
//				//image1f[tmpIndex] = (float)u1[tmpIndex];
//				//image2f[tmpIndex] = (float)u2[tmpIndex];
//
//				/*	if (nrhs > 6)
//				{
//				v1[tmpIndex] = (float)inputV[tmpIndex];
//				v2[tmpIndex] = (float)inputV[nPx + tmpIndex];
//				}
//				else
//				{
//				v1[tmpIndex] = 0.0f;
//				v2[tmpIndex] = 0.0f;
//				}
//				*/
//				/*v1[tmpIndex] = (float)inputV[tmpIndex];
//				v2[tmpIndex] = (float)inputV[nPx + tmpIndex];*/
//
//				/*		v1[tmpIndex] = (float)inputV(tmpIndex);
//				v2[tmpIndex] = (float)inputV(nPx + tmpIndex);*/
//				v1[tmpIndex] = _inputV(i, j, k, 0);
//				v2[tmpIndex] = _inputV(i, j, k, 1);
//				v3[tmpIndex] = _inputV(i, j, k, 2);
//				//_inputY1(i, j, 0, 0);
//				//v1[tmpIndex] = 0.0f;
//				//v2[tmpIndex] = 0.0f;
//
//
//				//cout << "v1[tmpIndex]: " << v1[tmpIndex] << " " << v2[tmpIndex] << " " << v3[tmpIndex] << endl;
//				Kty1[tmpIndex] = 0.0f;
//				Kty2[tmpIndex] = 0.0f;
//				Kty3[tmpIndex] = 0.0f;
//
//				//if (nrhs > 7)
//				//{
//				//	y11[tmpIndex] = (float)inputY[0 * nPx + tmpIndex];
//				//	y12[tmpIndex] = (float)inputY[1 * nPx + tmpIndex];
//				//	y21[tmpIndex] = (float)inputY[2 * nPx + tmpIndex];
//				//	y22[tmpIndex] = (float)inputY[3 * nPx + tmpIndex];
//				//             y3[tmpIndex] = (float)inputY[4 * nPx + tmpIndex];
//				//             
//				//             //for gradient constancy
//				//             y4[tmpIndex] = (float)inputY[5 * nPx + tmpIndex];
//				//             y5[tmpIndex] = (float)inputY[6 * nPx + tmpIndex];
//				//}
//				//else
//				//{
//				//	y11[tmpIndex] = 0.0f;
//				//	y12[tmpIndex] = 0.0f;
//				//	y21[tmpIndex] = 0.0f;
//				//	y22[tmpIndex] = 0.0f;
//				//             y3[tmpIndex] = 0.0f;
//				//             
//				//             //for gradient constancy
//				//             y4[tmpIndex] = 0.0f;
//				//             y5[tmpIndex] = 0.0f;
//				//}
//
//				/*	y11[tmpIndex] = (float)_inputY(0 * nPx + tmpIndex);
//				y12[tmpIndex] = (float)_inputY(1 * nPx + tmpIndex);
//				y21[tmpIndex] = (float)_inputY(2 * nPx + tmpIndex);
//				y22[tmpIndex] = (float)_inputY(3 * nPx + tmpIndex);
//				y3[tmpIndex] = (float)inputY(4 * nPx + tmpIndex);
//				*/
//				y11[tmpIndex] = _inputY1(i, j, k, 0);
//				y12[tmpIndex] = _inputY1(i, j, k, 1);
//				y13[tmpIndex] = _inputY1(i, j, k, 2);
//
//				y21[tmpIndex] = _inputY2(i, j, k, 0);
//				y22[tmpIndex] = _inputY2(i, j, k, 1);
//				y23[tmpIndex] = _inputY2(i, j, k, 2);
//
//				y31[tmpIndex] = _inputY3(i, j, k, 0);
//				y32[tmpIndex] = _inputY3(i, j, k, 1);
//				y33[tmpIndex] = _inputY3(i, j, k, 2);
//
//				y4[tmpIndex] = _inputY4(i, j, k);
//
//				//for gradient constancy
//				//	y4[tmpIndex] = (float)inputY(5 * nPx + tmpIndex);
//				//	y5[tmpIndex] = (float)inputY(6 * nPx + tmpIndex);
//
//
//
//				Kx11[tmpIndex] = 0.0f;
//				Kx12[tmpIndex] = 0.0f;
//				Kx13[tmpIndex] = 0.0f;
//
//				Kx21[tmpIndex] = 0.0f;
//				Kx22[tmpIndex] = 0.0f;
//				Kx23[tmpIndex] = 0.0f;
//
//				Kx31[tmpIndex] = 0.0f;
//				Kx32[tmpIndex] = 0.0f;
//				Kx33[tmpIndex] = 0.0f;
//
//				Kx4[tmpIndex] = 0.0f;
//
//				////for gradient constancy
//				//Kx4[tmpIndex] = 0.0f;
//				//Kx5[tmpIndex] = 0.0f;
//			}
//		}
//	}
//	cout << "Finish the initialization, open clock " << endl;
//	clock_t begin = clock();
//	//do k warpings
//	for (int idxwarp = 0; idxwarp < numberOfWarps; ++idxwarp)
//	{
//
//
//		// output??  u2, u21x,u21y, u22x,u22y
//		//		doWarp(  ut, ux, uy,uxt,uyt,uxx,uxy,uyx,uyy, sizeImage);
//		// output: ux, uy,  uxx,uxy,uyx,uyy
//		//  uxt,uyt : maybe for gradient consistency,
//		//ut,  
//
//		//doWarp(image1f, image2f, v1, v2, ut, ux, uy,uxt,uyt,uxx,uxy,uyx,uyy, sizeImage);
//		//doWarp(image1f, image2f, v1, v2, ut, ux, uy, sizeImage);
//
//
//		//////image test1(v1, sizeImage[0], sizeImage[1], 1, 1);
//		//////image test2(v2, sizeImage[0], sizeImage[1], 1, 1);
//
//		//real	m_aabb[] = { -584,584,-388,388,-10,10 };
//		//	real	m_aabb[] = { -584,584,-388,388,-10,10 };
//
//		//real	m_aabb[] = { -0.044,0.044,-0.045,0.045,-0.049,0.049 };
//		// 1. m_aabb issue/
//		// 2. try built-in cubic interpolation.
//		// 3. dt issue.
//		//	real	m_aabb[] = { -44,44,-45,45,-49,49 };
//		//real	m_aabb[] = { -100,100,-100,100,-100,100 };
//		real dt = 1;
//		//image _warpI2 = advection::advect(dim, _inputV, image2, -dt);
//		image _warpI2 = advection::advect(dim, _aabb, _inputV, image2, -dt);
//		//image _warpI2 = advection::advect(dim, _aabb, _inputV1, _inputV2, _inputV3, image2, -dt);
//		imagelist warpI2_xyz = _warpI2.get_gradient("xyz", 0);
//		////	image im_ux = warpI2_xy[0];
//		////	image im_uy = warpI2_xy[1];
//
//		//image ta = image2;
//		//real _max = image1.max();
//		//real _min= image1.min();
//		//image _warpI2 = image2.get_warp(_inputV, 1, 2, 0);//.normalize(_min, _max);
//		//image _warpI2 = image2.get_warp(_inputV, 1, 2, 0);//.normalize(_min, _max);
//
//		//image _warpI2 = image2.get_warp(_inputV, 1, 2, 0);//.normalize(_min, _max);
//		//imagelist warpI2_xyz = _warpI2.get_gradient("xyz", 0);
//
//
//		//#pragma omp parallel for
//		//for (int k = 0; k < dim[2]; ++k)
//		//{
//		//	for (int j = 0; j < dim[1]; ++j)
//		//	{
//		//		for (int i = 0; i < dim[0]; ++i)
//		//		{
//
//		#pragma omp parallel for
//		for (int k = ROI[4]; k < dim[2] - ROI[5]; ++k)
//		{
//			for (int j = ROI[2]; j < dim[1] - ROI[3]; ++j)
//			{
//				for (int i = ROI[0]; i < dim[0] - ROI[1]; ++i)
//				{
//					int tmpIndex = index3DtoLinear(dim, i, j, k);
//					//ux[tmpIndex] = im_ux(tmpIndex);
//					//uy[tmpIndex] = im_uy(tmpIndex);
//					////cout << ux[tmpIndex] << " " << uy[tmpIndex] << endl;
//					//ut[tmpIndex] = warpIm2(tmpIndex) - image1(tmpIndex) - ux[tmpIndex]*v1[tmpIndex] - uy[tmpIndex]*v2[tmpIndex];
//					//uXYZ(i,j,0,0)= warpI2_xy[0](i, j);
//					//uXYZ(i, j, 0, 1) = warpI2_xy[1](i, j);
//					ux[tmpIndex] = warpI2_xyz[0](i, j, k);
//					uy[tmpIndex] = warpI2_xyz[1](i, j, k);
//					uz[tmpIndex] = warpI2_xyz[2](i, j, k);
//					//cout << ux[tmpIndex] << " " << uy[tmpIndex] << " " << uz[tmpIndex] << endl;
//					//ut[tmpIndex] = _warpI2(i, j) - image1(i,j) - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex];
//					ut[tmpIndex] = _warpI2(i, j, k) - image1(i, j, k) - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex] - uz[tmpIndex] * v3[tmpIndex];
//					//ut[tmpIndex] = _warpI2(i, j) - image1(i, j) - uXYZ(i, j, 0, 0) * v1[tmpIndex] - uXYZ(i, j, 0, 1) * v2[tmpIndex];
//
//				}
//			}
//		}
//
//		//cout << " in the L1TVOpticalFlowNonlinear  2 " << endl;
//
//
//
//
//
//		#pragma omp parallel for
//		for (int i = 0; i < nPx; ++i)
//		{
//			//adjust step sizes
//			// GM Testing 1
//			tau1[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(ux[i]);
//			tau2[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(uy[i]);
//			tau3[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(uz[i]);
//			//tau1[i] = 4.0f / std::min(stepsize[1], stepsize[2]) + myAbs(uXYZ(i));
//			//tau2[i] = 4.0f / std::min(stepsize[1], stepsize[2]) + myAbs(uXYZ(i + nPx));
//			/*    if (gradientConstancy>0)
//			{
//			tau1[i] += std::abs(uxx[i]) + std::abs(uyx[i]);
//			tau2[i] += std::abs(uxy[i]) + std::abs(uyy[i]);
//			}*/
//			//sigma2[i] = std::abs(ux[i]) + std::abs(uXYZ(i + nPx));
//			sigma2[i] = std::abs(ux[i]) + std::abs(uy[i]) + std::abs(uz[i]);
//			//sigma2[i] = std::abs(uXYZ(i)) + std::abs(uXYZ(i+ nPx));
//
//			//for gradient constancy
//			/*         sigma3[i] = std::abs(uxx[i]) + std::abs(uxy[i]);
//			sigma4[i] = std::abs(uyx[i]) + std::abs(uyy[i]);
//			*/
//			tau1[i] = 1.0f / tau1[i];
//			tau2[i] = 1.0f / tau2[i];
//			tau3[i] = 1.0f / tau3[i];
//			sigma2[i] = 1.0f / sigma2[i];
//
//			//for gradient constancy
//			//sigma3[i] = 1.0f / sigma3[i];
//			//sigma4[i] = 1.0f / sigma4[i];
//		}
//
//
//		int iterations = 0;
//		float err = 1.0f;
//
//		while (err > tol && iterations <= maxIterations)
//		{
//			++iterations;
//
//			if (iterations % 50 == 0)
//			{
//				p = 0.0f;
//				d = 0.0f;
//			}
//
//			//cout << " in the L1TVOpticalFlowNonlinear  3  //primal step" << endl;
//
//			//primal step
//			//#pragma omp parallel for
//			//for (int k = 0; k < dim[2]; ++k)
//			//{
//			//	for (int j = 0; j < dim[1]; ++j)
//			//	{
//			//		for (int i = 0; i < dim[0]; ++i)
//			//		{
//			#pragma omp parallel for
//			for (int k = ROI[4]; k < dim[2] - ROI[5]; ++k)
//			{
//				for (int j = ROI[2]; j < dim[1] - ROI[3]; ++j)
//				{
//					for (int i = ROI[0]; i < dim[0] - ROI[1]; ++i)
//					{
//						int tmpIndex = index3DtoLinear(dim, i, j, k);
//
//						/*imagelist Y11grad = _inputY1.get_shared_channel(0).get_gradient("xy", -1);
//						imagelist Y12grad = _inputY1.get_shared_channel(1).get_gradient("xy", -1);
//						imagelist Y21grad = _inputY2.get_shared_channel(0).get_gradient("xy", -1);
//						imagelist Y22grad = _inputY2.get_shared_channel(1).get_gradient("xy", -1);*/
//						// we have some issues here.
//						/*#pragma omp parallel for
//						for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
//						{*/
//						if (iterations % 50 == 0)
//						{
//							v1Old[tmpIndex] = v1[tmpIndex];
//							v2Old[tmpIndex] = v2[tmpIndex];
//							v3Old[tmpIndex] = v3[tmpIndex];
//
//							Kty1Old[tmpIndex] = Kty1[tmpIndex];
//							Kty2Old[tmpIndex] = Kty2[tmpIndex];
//							Kty3Old[tmpIndex] = Kty3[tmpIndex];
//						}
//
//
//						//transpose equals -div  // m is mb  not choose 1
//						/*	Kty1[tmpIndex] = -stepsizeD[1] * Y11grad[0](tmpIndex)
//						- stepsizeD[2] * Y12grad[1](tmpIndex) + uXYZ(tmpIndex) * _inputY3(tmpIndex);
//
//						Kty2[tmpIndex] = -stepsizeD[1] * Y21grad[0](tmpIndex)
//						- stepsizeD[2] * Y22grad[1](tmpIndex) + uXYZ(tmpIndex + nPx)  * _inputY3(tmpIndex);*/
//
//						////transpose equals -div  // m is mb  not choose 1
//						//Kty1[tmpIndex] = -stepsizeD[1] * dxm(y11, dim, tableI[tmpIndex], tableJ[tmpIndex])
//						//	- stepsizeD[2] * dym(y12, dim, tableI[tmpIndex], tableJ[tmpIndex]) + ux[tmpIndex] * y3[tmpIndex];
//
//						//Kty2[tmpIndex] = -stepsizeD[1] * dxm(y21, dim, tableI[tmpIndex], tableJ[tmpIndex])
//						//	- stepsizeD[2] * dym(y22, dim, tableI[tmpIndex], tableJ[tmpIndex]) + uy[tmpIndex] * y3[tmpIndex];
//
//						//float dxm3( float *data,  int *sizeMat, int i, int j, int k)
//						//Kx11[tmpIndex] = stepsizeD[1] * dxp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//
//						//transpose equals -div  // m is mb  not choose 1
//						Kty1[tmpIndex] = -stepsizeD[1] * dxm3(y11, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
//							stepsizeD[2] * dym3(y12, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
//							stepsizeD[3] * dzm3(y13, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
//							ux[tmpIndex] * y4[tmpIndex];
//
//						Kty2[tmpIndex] = -stepsizeD[1] * dxm3(y21, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
//							stepsizeD[2] * dym3(y22, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
//							stepsizeD[3] * dzm3(y23, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
//							uy[tmpIndex] * y4[tmpIndex];
//
//						Kty3[tmpIndex] = -stepsizeD[1] * dxm3(y31, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
//							stepsizeD[2] * dym3(y32, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
//							stepsizeD[3] * dzm3(y33, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
//							uz[tmpIndex] * y4[tmpIndex];
//
//						//data[index3DtoLinear(dim, i - 1, j, k)]
//						//if (gradientConstancy>0)  // add a prior... means 
//						//{
//						//    Kty1[i] += uxx[i] * y4[i] + uyx[i] * y5[i];
//						//    Kty2[i] += uxy[i] * y4[i] + uyy[i] * y5[i];
//						//}
//						//cout << tau1[tmpIndex] << " " << tau2[tmpIndex] << " " << tau3[tmpIndex] << " " << endl;
//
//						//cout << Kty1[tmpIndex] << " " << Kty2[tmpIndex] << " " << Kty3[tmpIndex] << " " << endl;
//						v1[tmpIndex] -= tau1[tmpIndex] * Kty1[tmpIndex];
//						v2[tmpIndex] -= tau2[tmpIndex] * Kty2[tmpIndex];
//						v3[tmpIndex] -= tau3[tmpIndex] * Kty3[tmpIndex];
//						_inputV(i, j, k, 0) = v1[tmpIndex];// -= tau1[i] * Kty1[i];
//						_inputV(i, j, k, 1) = v2[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];
//						_inputV(i, j, k, 2) = v3[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];
//														   //if(v3[tmpIndex]>0.0f || v3[tmpIndex]<0.0f)
//														   //	cout << v1[tmpIndex] << " " << v2[tmpIndex] << " " << v3[tmpIndex] << " " << endl;
//														   //_inputV(tmpIndex) = v1[tmpIndex];// -= tau1[i] * Kty1[i];
//														   //_inputV(tmpIndex + nPx) = v2[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];
//
//
//														   //_inputV(tmpIndex) -= tau1[tmpIndex] * Kty1[tmpIndex];
//														   //_inputV(tmpIndex+ nPx) -= tau2[tmpIndex] * Kty2[tmpIndex];
//
//
//						if (iterations % 50 == 0)
//						{
//							//residuals
//							p += std::abs((v1Old[tmpIndex] - v1[tmpIndex]) / tau1[tmpIndex] - Kty1Old[tmpIndex] + Kty1[tmpIndex])
//								+ std::abs((v2Old[tmpIndex] - v2[tmpIndex]) / tau2[tmpIndex] - Kty2Old[tmpIndex] + Kty2[tmpIndex])
//								+ std::abs((v3Old[tmpIndex] - v3[tmpIndex]) / tau3[tmpIndex] - Kty3Old[tmpIndex] + Kty3[tmpIndex]);
//
//						}
//					}
//				}
//			}
//			//}		//cout << " in the L1TVOpticalFlowNonlinear  4  //dual step" << endl;
//
//
//			//dual step
//
//			//imagelist V1grad=_inputV.get_shared_channel(0).get_gradient("xy", 1);
//			//imagelist V2grad = _inputV.get_shared_channel(1).get_gradient("xy", 1);
//			/*imagelist V1grad = _inputV.get_shared_channel(0).get_gradient("xy", 1);
//			imagelist V2grad = _inputV.get_shared_channel(1).get_gradient("xy", 1);*/
//
//			//imagelist Vgrad = _inputV.get_gradient("xy", 1);
//			//image im1 = _inputV
//
//			#pragma omp parallel for reduction(+:d)
//			for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
//			{
//				float y11Tilde, y12Tilde, y13Tilde, y21Tilde, y22Tilde, y23Tilde, y31Tilde, y32Tilde, y33Tilde;
//
//				if (iterations % 50 == 0)
//				{
//
//
//					/*Y1old(tmpIndex) = _inputY1(tmpIndex);
//					Y1old(tmpIndex+ nPx) = _inputY1(tmpIndex+ nPx);
//
//					Y2old(tmpIndex) = _inputY2(tmpIndex);
//					Y2old(tmpIndex + nPx) = _inputY2(tmpIndex + nPx);
//					Y3old(tmpIndex) = _inputY3(tmpIndex);*/
//					y11Old[tmpIndex] = y11[tmpIndex];
//					y12Old[tmpIndex] = y12[tmpIndex];
//					y13Old[tmpIndex] = y13[tmpIndex];
//					y21Old[tmpIndex] = y21[tmpIndex];
//					y22Old[tmpIndex] = y22[tmpIndex];
//					y23Old[tmpIndex] = y23[tmpIndex];
//					y31Old[tmpIndex] = y31[tmpIndex];
//					y32Old[tmpIndex] = y32[tmpIndex];
//					y33Old[tmpIndex] = y33[tmpIndex];
//
//					y4Old[tmpIndex] = y4[tmpIndex];
//				}
//
//				/*KX1old(tmpIndex) = KX1(tmpIndex);
//				KX1old(tmpIndex+ nPx) = KX1(tmpIndex+ nPx);
//				KX2old(tmpIndex) = KX2(tmpIndex);
//				KX2old(tmpIndex + nPx) = KX2(tmpIndex + nPx);
//				KX3old(tmpIndex) = KX3(tmpIndex);*/
//				Kx11Old[tmpIndex] = Kx11[tmpIndex];
//				Kx12Old[tmpIndex] = Kx12[tmpIndex];
//				Kx13Old[tmpIndex] = Kx13[tmpIndex];
//				Kx21Old[tmpIndex] = Kx21[tmpIndex];
//				Kx22Old[tmpIndex] = Kx22[tmpIndex];
//				Kx23Old[tmpIndex] = Kx23[tmpIndex];
//				Kx31Old[tmpIndex] = Kx31[tmpIndex];
//				Kx32Old[tmpIndex] = Kx32[tmpIndex];
//				Kx33Old[tmpIndex] = Kx33[tmpIndex];
//
//				Kx4Old[tmpIndex] = Kx4[tmpIndex];
//				/*
//				KX1(tmpIndex) = stepsizeD[1] * V1grad[0](tmpIndex);
//				KX1(tmpIndex + nPx) = stepsizeD[2] * V1grad[1](tmpIndex);
//				KX2(tmpIndex) = stepsizeD[1] * V2grad[0](tmpIndex);
//				KX2(tmpIndex + nPx) = stepsizeD[2] * V2grad[1](tmpIndex);*/
//				/*Kx11[tmpIndex] = stepsizeD[1] * V1grad[0](tmpIndex);
//				Kx12[tmpIndex] = stepsizeD[2] * V1grad[1](tmpIndex);
//				Kx21[tmpIndex] = stepsizeD[1] * V2grad[0](tmpIndex);
//				Kx22[tmpIndex] = stepsizeD[2] * V2grad[1](tmpIndex);*/
//				/*	cout << Kx11[tmpIndex] << " " << Kx12[tmpIndex] << " " << Kx13[tmpIndex] << " " << Kx21[tmpIndex] << " " << Kx22[tmpIndex] << " " <<
//				Kx23[tmpIndex] << " " << endl;*/
//				Kx11[tmpIndex] = stepsizeD[1] * dxp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//				Kx12[tmpIndex] = stepsizeD[2] * dyp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//				Kx13[tmpIndex] = stepsizeD[3] * dzp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//
//				Kx21[tmpIndex] = stepsizeD[1] * dxp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//				Kx22[tmpIndex] = stepsizeD[2] * dyp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//				Kx23[tmpIndex] = stepsizeD[3] * dzp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//
//				Kx31[tmpIndex] = stepsizeD[1] * dxp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//				Kx32[tmpIndex] = stepsizeD[2] * dyp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//				Kx33[tmpIndex] = stepsizeD[3] * dzp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
//
//				Kx4[tmpIndex] = ux[tmpIndex] * v1[tmpIndex] + uy[tmpIndex] * v2[tmpIndex] + uz[tmpIndex] * v3[tmpIndex];  // warpped_u2=(ux, uy)
//																														  /*cout << Kx11[tmpIndex] << " " << Kx12[tmpIndex] << " " << Kx13[tmpIndex] << " " << Kx21[tmpIndex] << " " << Kx22[tmpIndex] << " " <<
//																														  Kx23[tmpIndex] << " " << endl;*/
//																														  //KX3(tmpIndex) = uXYZ(tmpIndex) * v1[tmpIndex] + uXYZ(tmpIndex + nPx) * v2[tmpIndex];  // warpped_u2=(ux, uy)
//																														  //Kx3[tmpIndex] = uXYZ(tmpIndex) * v1[tmpIndex] + uXYZ(tmpIndex+ nPx) * v2[tmpIndex];  // warpped_u2=(ux, uy)
//
//																														  /*Kx11[tmpIndex] = stepsizeD[1] * dxp(_inputV.get_channel(0), dim, tableI[tmpIndex], tableJ[tmpIndex]);
//																														  Kx12[tmpIndex] = stepsizeD[2] * dyp(_inputV.get_channel(0), dim, tableI[tmpIndex], tableJ[tmpIndex]);
//																														  Kx21[tmpIndex] = stepsizeD[1] * dxp(_inputV.get_channel(1), dim, tableI[tmpIndex], tableJ[tmpIndex]);
//																														  Kx22[tmpIndex] = stepsizeD[2] * dyp(_inputV.get_channel(1), dim, tableI[tmpIndex], tableJ[tmpIndex]);*/
//
//																														  //Kx3[tmpIndex] = ux[tmpIndex] * _inputV(tmpIndex) + uy[tmpIndex] * _inputV(tmpIndex+ nPx);  // warpped_u2=(ux, uy)
//
//																														  /*if (gradientConstancy>0)
//																														  {
//																														  Kx4Old[tmpIndex] = Kx4[tmpIndex];
//																														  Kx5Old[tmpIndex] = Kx5[tmpIndex];
//
//																														  Kx4[tmpIndex] = uxx[tmpIndex] * v1[tmpIndex] + uxy[tmpIndex] * v2[tmpIndex];
//																														  Kx5[tmpIndex] = uyx[tmpIndex] * v1[tmpIndex] + uyy[tmpIndex] * v2[tmpIndex];
//
//																														  y4Old[tmpIndex] = y4[tmpIndex];
//																														  y5Old[tmpIndex] = y5[tmpIndex];
//
//																														  y4[tmpIndex] = std::max(-1.0f,std::min(1.0f, y4[tmpIndex] + sigma3[tmpIndex]*(Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old[tmpIndex] + uxt[tmpIndex])));
//																														  y5[tmpIndex] = std::max(-1.0f,std::min(1.0f, y5[tmpIndex] + sigma4[tmpIndex]*(Kx5[tmpIndex] + Kx5[tmpIndex] - Kx5Old[tmpIndex] + uyt[tmpIndex])));
//																														  }*/
//
//				if (typeNorm == 4) // Huber
//				{
//					/*y11Tilde = (_inputY1(tmpIndex) + sigma1*(KX1(tmpIndex) + KX1(tmpIndex) - KX1old(tmpIndex))) * huberFactor;
//					y12Tilde = (_inputY1(tmpIndex + nPx) + sigma1*(KX1(tmpIndex + nPx) + KX1(tmpIndex + nPx) - KX1old(tmpIndex + nPx))) * huberFactor;
//					y21Tilde = (_inputY2(tmpIndex)+ sigma1*(KX2(tmpIndex) + KX2(tmpIndex) - KX2old(tmpIndex))) * huberFactor;
//					y22Tilde = (_inputY2(tmpIndex + nPx) + sigma1*(KX2(tmpIndex + nPx) + KX2(tmpIndex + nPx) - KX2old(tmpIndex + nPx))) * huberFactor;*/
//					y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex])) * huberFactor;
//					y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex])) * huberFactor;
//					y13Tilde = (y13[tmpIndex] + sigma1*(Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old[tmpIndex])) * huberFactor;
//
//					y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex])) * huberFactor;
//					y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex])) * huberFactor;
//					y23Tilde = (y23[tmpIndex] + sigma1*(Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old[tmpIndex])) * huberFactor;
//
//					y31Tilde = (y31[tmpIndex] + sigma1*(Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old[tmpIndex])) * huberFactor;
//					y32Tilde = (y32[tmpIndex] + sigma1*(Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old[tmpIndex])) * huberFactor;
//					y33Tilde = (y33[tmpIndex] + sigma1*(Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old[tmpIndex])) * huberFactor;
//				}
//				else
//				{
//					/*y11Tilde = _inputY1(tmpIndex) + sigma1*(KX1(tmpIndex) + KX1(tmpIndex) - KX1old(tmpIndex));
//					y12Tilde = _inputY1(tmpIndex + nPx) + sigma1*(KX1(tmpIndex + nPx) + KX1(tmpIndex + nPx) - KX1old(tmpIndex + nPx));
//					y21Tilde = _inputY2(tmpIndex) + sigma1*(KX2(tmpIndex) + KX2(tmpIndex) - KX2old(tmpIndex));
//					y22Tilde = _inputY2(tmpIndex + nPx) + sigma1*(KX2(tmpIndex + nPx) + KX2(tmpIndex + nPx) - KX2old(tmpIndex + nPx)) ;*/
//					// ATV or ITV?
//
//					y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex]));
//					y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex]));
//					y13Tilde = (y13[tmpIndex] + sigma1*(Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old[tmpIndex]));
//
//					y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex]));
//					y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex]));
//					y23Tilde = (y23[tmpIndex] + sigma1*(Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old[tmpIndex]));
//
//					y31Tilde = (y31[tmpIndex] + sigma1*(Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old[tmpIndex]));
//					y32Tilde = (y32[tmpIndex] + sigma1*(Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old[tmpIndex]));
//					y33Tilde = (y33[tmpIndex] + sigma1*(Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old[tmpIndex]));
//				}
//
//				float divisor1 = std::max(1.0f, std::sqrt(y11Tilde*y11Tilde + y12Tilde*y12Tilde + y13Tilde*y13Tilde) / lambda);
//				float divisor2 = std::max(1.0f, std::sqrt(y21Tilde*y21Tilde + y22Tilde*y22Tilde + y23Tilde*y23Tilde) / lambda);
//				float divisor3 = std::max(1.0f, std::sqrt(y31Tilde*y31Tilde + y32Tilde*y32Tilde + y33Tilde*y33Tilde) / lambda);
//
//
//				//cout << "divisor2: "<< divisor1 << " " << divisor2 << " " << divisor3 << endl;
//				/*_inputY1(tmpIndex) = y11Tilde / divisor1;
//				_inputY1(tmpIndex+ nPx) = y12Tilde / divisor1;
//				_inputY2(tmpIndex) = y21Tilde / divisor2;
//				_inputY2(tmpIndex+ nPx) = y22Tilde / divisor2;*/
//				y11[tmpIndex] = y11Tilde / divisor1;
//				y12[tmpIndex] = y12Tilde / divisor1;
//				y13[tmpIndex] = y13Tilde / divisor1;
//
//				y21[tmpIndex] = y21Tilde / divisor2;
//				y22[tmpIndex] = y22Tilde / divisor2;
//				y23[tmpIndex] = y23Tilde / divisor2;
//
//				y31[tmpIndex] = y31Tilde / divisor3;
//				y32[tmpIndex] = y32Tilde / divisor3;
//				y33[tmpIndex] = y33Tilde / divisor3;
//
//
//
//				//_inputY3(tmpIndex) = std::max(-1.0f, std::min(1.0f, _inputY3(tmpIndex) + sigma2[tmpIndex] * (KX3(tmpIndex) + KX3(tmpIndex) - KX3old(tmpIndex) + ut[tmpIndex])));
//
//				// where to get the u_t, 
//				//y3[tmpIndex] = std::max(-1.0f, std::min(1.0f, y3[tmpIndex] + sigma2[tmpIndex] * (Kx3[tmpIndex] + Kx3[tmpIndex] - Kx3Old[tmpIndex] + ut[tmpIndex])));
//				y4[tmpIndex] = std::max(-1.0f, std::min(1.0f, y4[tmpIndex] + sigma2[tmpIndex] * (Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old[tmpIndex] + ut[tmpIndex])));
//
//				if (iterations % 50 == 0)
//				{
//
//					//d += std::abs((Y1old(tmpIndex) - _inputY1(tmpIndex)) / sigma1 - KX1old(tmpIndex) + KX1(tmpIndex)) +
//					//	std::abs((Y1old(tmpIndex + nPx) - _inputY1(tmpIndex + nPx)) / sigma1 - KX1old(tmpIndex + nPx) + KX1(tmpIndex + nPx)) +
//					//	std::abs((Y2old(tmpIndex) - _inputY2(tmpIndex)) / sigma1 - KX2old(tmpIndex) + KX2(tmpIndex)) +
//					//	std::abs((Y2old(tmpIndex + nPx) - _inputY2(tmpIndex + nPx)) / sigma1 - KX2old(tmpIndex + nPx) + KX2(tmpIndex + nPx)) +
//					//	std::abs((Y3old(tmpIndex) - _inputY3(tmpIndex)) / sigma2[tmpIndex] - KX3old(tmpIndex) + KX3(tmpIndex));
//
//					d += std::abs((y11Old[tmpIndex] - y11[tmpIndex]) / sigma1 - Kx11Old[tmpIndex] + Kx11[tmpIndex]) +
//						std::abs((y12Old[tmpIndex] - y12[tmpIndex]) / sigma1 - Kx12Old[tmpIndex] + Kx12[tmpIndex]) +
//						std::abs((y13Old[tmpIndex] - y13[tmpIndex]) / sigma1 - Kx13Old[tmpIndex] + Kx13[tmpIndex]) +
//
//						std::abs((y21Old[tmpIndex] - y21[tmpIndex]) / sigma1 - Kx21Old[tmpIndex] + Kx21[tmpIndex]) +
//						std::abs((y22Old[tmpIndex] - y22[tmpIndex]) / sigma1 - Kx22Old[tmpIndex] + Kx22[tmpIndex]) +
//						std::abs((y23Old[tmpIndex] - y23[tmpIndex]) / sigma1 - Kx23Old[tmpIndex] + Kx23[tmpIndex]) +
//
//						std::abs((y31Old[tmpIndex] - y31[tmpIndex]) / sigma1 - Kx31Old[tmpIndex] + Kx31[tmpIndex]) +
//						std::abs((y32Old[tmpIndex] - y32[tmpIndex]) / sigma1 - Kx32Old[tmpIndex] + Kx32[tmpIndex]) +
//						std::abs((y33Old[tmpIndex] - y33[tmpIndex]) / sigma1 - Kx33Old[tmpIndex] + Kx33[tmpIndex]) +
//
//						std::abs((y4Old[tmpIndex] - y4[tmpIndex]) / sigma2[tmpIndex] - Kx4Old[tmpIndex] + Kx4[tmpIndex]);
//
//					/*if (gradientConstancy>0)
//					{
//					d += std::abs((y4Old[tmpIndex] - y4[tmpIndex]) / sigma3[tmpIndex] - Kx4Old[tmpIndex] + Kx4[tmpIndex]) +
//					std::abs((y5Old[tmpIndex] - y5[tmpIndex]) / sigma4[tmpIndex] - Kx5Old[tmpIndex] + Kx5[tmpIndex]);
//					}*/
//				}
//			}
//
//			if (iterations % 50 == 0)
//			{
//				err = (d*d + p*p) / (float)nPx;
//			}
//
//			if (iterations % 1000 == 0)
//			{
//
//				cout << "Iteration: " << iterations << " Residual " << err << endl;
//
//				//mexPrintf("Iteration %d,Residual %e\n", iterations, err);
//				//mexEvalString("drawnow;");
//			}
//		}
//	}
//	//cout << " in the L1TVOpticalFlowNonlinear 5 " << endl;
//
//	//write output
//	//#pragma omp parallel for
//	//for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
//	//{
//
//	// end of warps 
//
//	//#pragma omp parallel for // collapse(3)
//	//for (int k = 0; k < dim[2]; ++k)
//	//{
//	//	for (int j = 0; j < dim[1]; ++j)
//	//	{
//	//		for (int i = 0; i < dim[0]; ++i)
//	//		{
//
//	#pragma omp parallel for
//	for (int k = ROI[4]; k < dim[2] - ROI[5]; ++k)
//	{
//		for (int j = ROI[2]; j < dim[1] - ROI[3]; ++j)
//		{
//			for (int i = ROI[0]; i < dim[0] - ROI[1]; ++i)
//			{
//				int tmpIndex = index3DtoLinear(dim, i, j, k);
//				/*
//				YOut[tmpIndex + 0 * nPx] = (double)y11[tmpIndex];
//				YOut[tmpIndex + 1 * nPx] = (double)y12[tmpIndex];
//				YOut[tmpIndex + 2 * nPx] = (double)y21[tmpIndex];
//				YOut[tmpIndex + 3 * nPx] = (double)y22[tmpIndex];
//				YOut[tmpIndex + 4 * nPx] = (double)y3[tmpIndex];
//
//				YOut[tmpIndex + 5 * nPx] = (double)y4[tmpIndex];
//				YOut[tmpIndex + 6 * nPx] = (double)y5[tmpIndex];
//
//				Outv1[tmpIndex] = (double)v1[tmpIndex];
//				Outv2[tmpIndex] = (double)v2[tmpIndex];*/
//
//				_inputY1(i, j, k, 0) = (float)y11[tmpIndex];
//				_inputY1(i, j, k, 1) = (float)y12[tmpIndex];
//				_inputY1(i, j, k, 2) = (float)y13[tmpIndex];
//
//				_inputY2(i, j, k, 0) = (float)y21[tmpIndex];
//				_inputY2(i, j, k, 1) = (float)y22[tmpIndex];
//				_inputY2(i, j, k, 2) = (float)y23[tmpIndex];
//
//				_inputY3(i, j, k, 0) = (float)y31[tmpIndex];
//				_inputY3(i, j, k, 1) = (float)y32[tmpIndex];
//				_inputY3(i, j, k, 2) = (float)y33[tmpIndex];
//
//				_inputY4(i, j, k, 0) = (float)y4[tmpIndex];
//
//				//	_inputY(tmpIndex + 5 * nPx) = (double)y4[tmpIndex];
//				//	_inputY(tmpIndex + 6 * nPx) = (double)y5[tmpIndex];
//
//				//	_inputV[tmpIndex] += (double)v1[tmpIndex];
//				//	_inputV[tmpIndex +  nPx] += (double)v2[tmpIndex];
//				_inputV(i, j, k, 0) = (float)v1[tmpIndex];
//				_inputV(i, j, k, 1) = (float)v2[tmpIndex];
//				_inputV(i, j, k, 2) = (float)v3[tmpIndex];
//				//	cout << _inputV(i, j, k, 0) << " " << _inputV(i, j, k, 1) << " " << _inputV(i, j, k, 2) << endl;
//
//				/*		Outv1[tmpIndex] = _inputV(i, j, k, 0);
//				Outv2[tmpIndex] = _inputV(i, j, k, 1);
//				Outv3[tmpIndex] = _inputV(i, j, k, 2);*/
//				//	_inputV(tmpIndex)
//				//cout << Outv1[tmpIndex] << " ";
//			}
//
//		}
//	}
//	cout << " in the L1TVOpticalFlowNonlinear  5 " << endl;
//
//	clock_t end = clock();
//	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	cout << "=========Total time is :" << elapsed_secs << " s ======" << endl;
//	int num1 = _inputV.width();
//	int num2 = _inputV.height();
//	int num3 = _inputV.depth();
//	//itoa(num1, snum1, 10);
//	//itoa(num2, snum2, 10);
//	//str = snum1+"h" + snum2+"-1";
//	//str+="h"+
//
//	char out_dir[1024];
//	char basename[1024];
//	char flow_field[1024];
//	//sprintf(out_dir, ".");
//	sprintf(basename, "outputf1f2Roi");
//	//sprintf(flow_field, "flow_field");
//
//
//	//char Mybuf[10];
//	//sprintf(Mybuf, "%d", num1);
//	//string Myb = Mybuf;
//
//	//char Mybuf2[10];
//	//sprintf(Mybuf2, "%d", num2);
//	//string Myb2 = Mybuf2;
//
//	//char Mybuf3[10];
//	//sprintf(Mybuf3, "%d", num3);
//	//string Myb3 = Mybuf3;
//
//	//string _output;// = pOutput;
//	//_output += "/";
//
//	//string Mywritefile1 = Myb + "-" + Myb2 + "-" + Myb3 + "-1.txt";
//	//string Mywritefile2 = Myb + "-" + Myb2 + "-" + Myb3 + "-2.txt";
//	//string Mywritefile3 = Myb + "-" + Myb2 + "-" + Myb3 + "-3.txt";
//
//	//string Mywritefile4 = Myb + "-" + Myb2 + "-" + Myb3 + "-volume.raw";
//	//string Mywritefile5 = "cubic/" + Myb + "-" + Myb2 + "-" + Myb3 + "-vol.hdr";
//	//string Mywritefile6 = "cubic/" + Myb + "-" + Myb2 + "-" + Myb3 + "-field";
//
//
//	//string file1 = "cubic/"+Myb + "-" + Myb2 + "-" + Myb3 + "backward1dt.hdr";
//	//string file2 = "cubic/" + Myb + "-" + Myb2 + "-" + Myb3 + "forward1dt.hdr";
//	//string file3 = "cubic/" + Myb + "-" + Myb2 + "-" + Myb3 + "vol1-sample.hdr";
//	//string file4 = "cubic/" + Myb + "-" + Myb2 + "-" + Myb3 + "vol2-sample.hdr";
//
//
//	//char* MyOutputfile1 = _strdup(Mywritefile1.c_str());
//	//char* MyOutputfile2 = _strdup(Mywritefile2.c_str());
//	//char* MyOutputfile3 = _strdup(Mywritefile3.c_str());
//	//char* MyOutputfile4 = _strdup(Mywritefile4.c_str());
//	//char* MyOutputfile5 = _strdup(Mywritefile5.c_str());
//	//char* MyOutputfile6 = _strdup(Mywritefile6.c_str());
//	//char* _file1 = _strdup(file1.c_str());
//	//char* _file2 = _strdup(file2.c_str());
//	//char* _file3 = _strdup(file3.c_str());
//	//char* _file4 = _strdup(file4.c_str());
//	//strdup(
//	//FILE *fp;
//	////fp = fopen("test1.txt", "wb");
//	//fp = fopen(MyOutputfile1, "wb");
//
//	///*if (fp == null)
//	//return;*/
//	////char x[10] = "ABCDEFGHIJ";
//	//fwrite(Outv1, sizeof(float), dim[0] * dim[1] * dim[2], fp);
//	//fclose(fp);
//
//
//	//FILE *fp2;
//	//fp2 = fopen(MyOutputfile2, "wb");
//	///*if (fp == null)
//	//return;*/
//	////char x[10] = "ABCDEFGHIJ";
//	//fwrite(Outv2, sizeof(float), dim[0] * dim[1] * dim[2], fp2);
//	//fclose(fp2);
//
//	//FILE *fp3;
//	//fp3 = fopen(MyOutputfile3, "wb");
//	///*if (fp == null)
//	//return;*/
//	////char x[10] = "ABCDEFGHIJ";
//	//fwrite(Outv3, sizeof(float), dim[0] * dim[1] * dim[2], fp3);
//	//fclose(fp3);
//
//
//
//
//	//	real	m_aabb[] = { -584,584,-388,388,-10,10 };
//
//	//real	m_aabb[] = { -0.044,0.044,-0.045,0.045,-0.049,0.049 };
//	//real	m_aabb[] = { -44,44,-45,45,-49,49 };
//
//	real dt = 1;
//	image _warpI2 = advection::advect(dim, _aabb, _inputV, image2, -dt);
//	//_warpI2.save_raw(MyOutputfile5);
//	//_warpI2.normalize(image1.normalize(0,40.0f);
//	//image2.normalize(0, 40.0f);)
//	//_warpI2.save_analyze(MyOutputfile5);
//	//str_format("%s/%s.%04d.vtr",out_dir,basename,iter).c_str()
//	// out . 
//	// base  output
//
//
//	_warpI2.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.hdr", basename, num1, num2, num3, _lambda, huber).c_str());
//
//	//image T_warpI2 = image2.get_warp(_inputV, 1, 2, 0);//.normalize(_min, _max);
//	//T_warpI2.save_analyze("built-inWarp_I2backward_absolute.hdr");
//	//imagelist warpI2_xyz = _warpI2.get_gradient("xyz", 0);
//
//	//image T_warpI3 = image1.get_warp(_inputV, 3, 2, 0);//.normalize(_min, _max);
//	//T_warpI3.save_analyze("built-inWarp_I1forward_absolute.hdr");
//	//float *testData = _warpI2;
//	//_warpI2.save_empty_cimg("testfile_warpI2.cimg", _warpI2.width(), _warpI2.height(), _warpI2.depth(), _warpI2.spectrum());
//	//_inputV.save_empty_cimg(MyOutputfile6, _inputV.width(), _inputV.height(), _inputV.depth(), _inputV.spectrum());
//	//FILE *fp4;
//	//fp4= fopen(MyOutputfile4, "wb");
//	//m_Data->m_Volumelist[iframes].save_analyze(MyOutputfile2);
//
//	///*if (fp == null)
//	//return;*/
//	////char x[10] = "ABCDEFGHIJ";
//	//fwrite(testData, sizeof(float), dim[0] * dim[1] * dim[2], fp4);
//	//fclose(fp4);
//	//Sleep(1000);
//	cout << "Testing the application end part" << endl;
//
//
//
//	//real	m_aabb[] = { -584,584,-388,388,-10,10 };
//	//dt = 5;
//	//image _warpI3 = advection::advect(dim, m_aabb, _inputV, image2, -dt);
//	//_warpI3.save_raw("Positive_5_Volume.cimg");
//
//	//m_Data->m_Volumelist[iframes].save_analyze(MyOutputfile2);
//
//	//dt = -5;
//	//image _warpI4 = advection::advect(dim, m_aabb, _inputV, image2, -dt);
//	//_warpI4.save_raw("Minues_5_Volume.cimg");
//
//	_inputV.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.flowfield.hdr", basename, num1, num2, num3, _lambda, huber).c_str());
//	//_inputV.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.flowfield.hdr", basename, num1, num2, num3, _lambda, huber).c_str());
//	//_inputV.save_empty_cimg(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.flowfield_test.hdr", basename, num1, num2, num3, _lambda, huber).c_str(), _inputV.width(), _inputV.height(), _inputV.depth(), _inputV.spectrum());
//	_inputV.get_channel(0).save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.flowV1.hdr", basename, num1, num2, num3, _lambda, huber).c_str());
//	_inputV.get_channel(1).save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.flowV2.hdr", basename, num1, num2, num3, _lambda, huber).c_str());
//	_inputV.get_channel(2).save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.flowV3.hdr", basename, num1, num2, num3, _lambda, huber).c_str());
//	//_inputV.save_empty_cimg(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.flowfield_v1.hdr", basename, num1, num2, num3, _lambda, huber).c_str(), _inputV.width(), _inputV.height(), _inputV.depth(), _inputV.spectrum());
//
//	dt = 1;
//	image _warpI5 = advection::advect(dim, _aabb, _inputV, image2, -dt);
//	//_warpI5.save_raw("forward_FromVolume1.cimg");
//	_warpI5.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.backward.hdr", basename, num1, num2, num3, _lambda, huber).c_str());
//
//	//_warpI5.normalize(image1.min(), image1.max());
//	//_warpI5.save_analyze(_file1);
//
//
//	dt = -1;
//	image _warpI6 = advection::advect(dim, _aabb, _inputV, image1, -dt);
//	_warpI6.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.forward.%5.2fs.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str());
//
//	/*dt = -1.01;
//	image _warpI8 = advection::advect(dim, _aabb, _inputV, image1, -dt);
//	_warpI8.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.forward1.04dt.%5.2fs.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str());*/
//
//
//	//dt = -1.02;
//	//image _warpI7 = advection::advect(dim, _aabb, _inputV, image1, -dt);
//	//_warpI7.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.forward1.05dt.%5.2fs.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str());
//
//
//	//dt = -1.06;
//	//image _warpI9 = advection::advect(dim, _aabb, _inputV, image1, -dt);
//	//_warpI9.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.forward1.06dt.%5.2fs.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str());
//
//	/*dt = -1.1;
//	image _warpI9 = advection::advect(dim, _aabb, _inputV, image1, -dt);
//	_warpI9.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.forward1.2dt.%5.2fs.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str());
//	*/
//	//_warpI6.normalize(image2.min(), image2.max());
//	//_warpI6.save_analyze(_file2);
//
//	/*dt = 10;
//	image _warpI7 = advection::advect(dim, _aabb, _inputV, image2, -dt);
//	_warpI7.save_analyze("cubic/Backward_10dt.hdr");*/
//
//	//dt = -1;
//	//image _warpI7 = advection::advect(dim, _aabb, _inputV, image1, -dt);
//	////_warpI5.save_raw("forward_FromVolume1.cimg");
//	//_warpI7.save_analyze("Forward_FromVolume1_1dt.hdr");
//
//
//	//dt = -5;
//	//image _warpI8 = advection::advect(dim, _aabb, _inputV, image1, -dt);
//	//_warpI8.save_analyze("Forward_FromVolume1_5dt.hdr");
//	//str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.forward1.05dt.%5.2fs.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str()
//	/*image1.save_analyze(_file3);
//	image2.save_analyze(_file4);*/
//	image1.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.sample1.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str());
//	image2.save_analyze(str_format("%s/%03d.%03d.%03d.s%3.3f.h%3.3f.sample2.hdr", basename, num1, num2, num3, _lambda, huber, elapsed_secs).c_str());
//	//_warpI6.save_raw("forward_FromVolume1_dt5.cimg");
//	//dt = 0.1;
//	//image _warpI5 = advection::advect(dim, m_aabb, _inputV, image2, -dt);
//	//_warpI5.save_raw("Minues_0.1_Volume.cimg");
//
//	//image warpField2(sizeImage[0], sizeImage[1], 1, 2);
//
//	////cimg_forXY(warpField, x, y) {
//	//for (int j = 0; j < sizeImage[1]; j++)
//	//{
//	//	for (int i = 0; i < sizeImage[0]; i++)
//	//	{
//	//		int tmpIndex = index2DtoLinear(sizeImage, i, j);
//	//		//const float
//	//		//u = x + x0 - dw2, v = y + y0 - dh2;
//
//	//		warpField2(i, j, 0) = -Outv1[tmpIndex];
//	//		warpField2(i, j, 1) = -Outv2[tmpIndex];
//	//	}
//	//}
//	//image ta = image2;
//	//real _max = image1.max();
//	//real _min= image1.min();
//	//image warpIm22 = ta.get_warp(warpField2, 1, 2, 0).normalize(_min, _max);
//
//
//	//image test2(Outv1, sizeImage[0], sizeImage[1], 1, 1);
//	//image test3(Outv2, sizeImage[0], sizeImage[1], 1, 1);
//	////test2 = -test2;
//	////test3 = -test3;
//	////cout << i2(0) << " " << i2(4) << " " << i2(8) << endl;
//	//image testwarpfield2 = test2.get_append(test3,'c');
//	////testwarpfield2.insert(test2);
//	////testwarpfield2.insert(test3);
//	////image testfiled = testwarpfield2.get_append();
//	//image tt = image2;
//	//image warpIm2t = tt.get_warp(testwarpfield2, 1, 2,1);
//
//	////image testwarpfield = image(v1).get_append(testymap, 'c');
//	/*
//	imagelist testresultlist(image2, image1, warpIm22);
//
//	Sleep(500);
//
//	testresultlist.display();
//	Sleep(500);
//	*/
//
//	/*	CImgList<float> teste(abs(image2 - warpIm2), abs(image2 - warpIm2));
//	teste.display();*/
//
//
//
//	delete[] tableI;
//	delete[] tableJ;
//	delete[] tableK;
//
//	delete[] image1f;
//	delete[] image2f;
//
//	delete[] sigma2;
//	//delete[] sigma3;
//	//delete[] sigma4;
//	delete[] tau1;
//	delete[] tau2;
//	delete[] tau3;
//
//	delete[] v1;
//	delete[] v2;
//	delete[] v3;
//	delete[] v1Old;
//	delete[] v2Old;
//	delete[] v3Old;
//
//	delete[] ux;
//	delete[] uy;
//	delete[] uz;
//	delete[] ut;
//
//	//delete[] uxx;
//	//delete[] uxy;
//	//delete[] uyx;
//	//delete[] uyy;
//
//	//delete[] uxt;
//	//delete[] uyt;
//
//	delete[] y11;
//	delete[] y12;
//	delete[] y13;
//	delete[] y21;
//	delete[] y22;
//	delete[] y23;
//	delete[] y31;
//	delete[] y32;
//	delete[] y33;
//	delete[] y4;
//	//delete[] y4;
//	//delete[] y5;
//
//	delete[] y11Old;
//	delete[] y12Old;
//	delete[] y13Old;
//	delete[] y21Old;
//	delete[] y22Old;
//	delete[] y23Old;
//	delete[] y31Old;
//	delete[] y32Old;
//	delete[] y33Old;
//	delete[] y4Old;
//	//delete[] y4Old;
//	//delete[] y5Old;
//
//	delete[] Kty1;
//	delete[] Kty2;
//	delete[] Kty3;
//
//	delete[] Kty1Old;
//	delete[] Kty2Old;
//	delete[] Kty3Old;
//
//	delete[] Kx11;
//	delete[] Kx12;
//	delete[] Kx13;
//	delete[] Kx21;
//	delete[] Kx22;
//	delete[] Kx23;
//	delete[] Kx31;
//	delete[] Kx32;
//	delete[] Kx33;
//	delete[] Kx4;
//	//delete[] Kx3;
//	//delete[] Kx4;
//	//delete[] Kx5;
//
//	delete[] Kx11Old;
//	delete[] Kx12Old;
//	delete[] Kx13Old;
//	delete[] Kx21Old;
//	delete[] Kx22Old;
//	delete[] Kx23Old;
//	delete[] Kx31Old;
//	delete[] Kx32Old;
//	delete[] Kx33Old;
//	delete[] Kx4Old;
//	//delete[] Kx4Old;
//	//delete[] Kx5Old;
//
//	//return ;
//}
//
//
//
//
//
//
////template< typename real, typename image >
//void L1TVOpticalFlowNonlinear_MultiScale3D( int *dim, const image &image1, const image &image2, image &_inputV,
//	image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4, real _tol, real _lambda, int _maxIterations, int _norm, int _numberOfWarps,
//	real eta, real sigma, int scales, real _huber)
//{
//
//
//	//m_aabb[0] = -dim[0] * 100;
//	//m_aabb[1] = dim[0]*100 ;
//	//m_aabb[2] = -dim[1] * 100;
//	//m_aabb[3] = dim[1]*100;
//	//m_aabb[4] = -dim[2] * 100;
//	//m_aabb[5] = dim[2]*100;
//	//if (levels > 1)
//	//{
//
//	//	std::cout << "V size: " << _inputV.width() << " " << _inputV.height() << std::endl;
//	//	image I1_coarse = blur_and_downsample(image1, eta, sigma);
//	//	image I2_coarse = blur_and_downsample(image2, eta, sigma);
//	//	image V_coarse = Y_blur_and_downsample(_inputV, eta, sigma);
//	//	// maybe need to seperate the y1, y2, y
//	//	//image Y_coarse= blur_and_downsample(_inputY, eta, sigma);
//	//	image Y_coarse = Y_blur_and_downsample(_inputY, eta, sigma);
//	//	image V_coarse_init = V_coarse;
//	//	int tdim[] = { I1_coarse.width(), I1_coarse.height() };
//	//	std::cout << "tdim size: " << tdim[0] <<" "<<tdim[1]<< std::endl;
//
//	//int tdim[] = { I1_coarse.width(), I1_coarse.height(), I1_coarse.depth() };
//	//	//optical_flow_horn_schunk_multiscale(tdim, aabb, I1_coarse, I2_coarse, uvw_coarse, opts, point_data_coarse, cell_data_coarse, level - 1);
//	//
//
//	//	 L1TVOpticalFlowNonlinear_MultiScale(tdim, I1_coarse, I2_coarse, V_coarse,
//	//										 Y_coarse, _tol, _lambda, _maxIterations, _norm, _numberOfWarps,eta,  sigma,  levels-1);
//
//
//
//	//	 V_coarse += (V_coarse - V_coarse_init).get_resize(dim[0], dim[1], 2,3, 3);
//
//	//}
//	//std::cout <<  "optical flow level: " << levels << std::endl;
//
//
//	//return L1TVOpticalFlowNonlinear(dim, image1, image2, _inputV, _inputY, _tol, _lambda, _maxIterations, _norm,
//	//								_numberOfWarps);
//
//
//	//image I1_coarse = blur_and_downsample(image1, eta, sigma);
//	//image I2_coarse = blur_and_downsample(image2, eta, sigma);
//	//image V_coarse = Y_blur_and_downsample(_inputV, eta, sigma);
//	image I1_blur = image1.get_blur(sigma);
//	image I2_blur = image2.get_blur(sigma);
//
//	/*image tempI1 = I1_blur;
//	image tempI2 = I2_blur;
//	image tempV = _inputV;
//	image tempY1 = _inputY1;
//	image tempY2 = _inputY2;
//	image tempY3 = _inputY3;*/
//
//	//image I3_blur = image3.get_blur(sigma);
//	//float **Im1s;
//	//float **Im2s;
//	int size = image1.width()*image1.height()*image1.depth();
//	//Im1s[0] = new float[size];
//	//Im2s[0] = new float[size];
//	imagelist Im1s;
//	imagelist Im2s;
//	imagelist inputVs;
//	//imagelist inputV2s;
//	//imagelist inputV3s;
//	imagelist inputY1s;
//	imagelist inputY2s;
//	imagelist inputY3s;
//	imagelist inputY4s;
//
//	//// load old version
//
//
//	////imagelist dims;
//	//int twidth = I1_blur.width();
//	//int theight = I1_blur.height();
//	//int tdepth = I1_blur.depth();
//	//float t_eta = 1.0f;
//
//	//int tdims[] = { twidth ,theight, tdepth };
//
//	//cout << tdims[0] << " " << tdims[1] << " " << tdims[2] << endl;
//	//for (int i = 0; i < scales; i++)
//	//{
//
//	//	Im1s.insert(I1_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//	//	//Im1s.insert()
//	//	//image timg = tempI1.get_resize(tempI1.width()*eta, tempI1.height()*eta, tempI1.depth()*eta, 1, 5);
//	//	//tempI1 = timg;
//	//	Im2s.insert(I2_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//	//	//Im2s.insert(tempI2);
//	//	//image timg2 = tempI2.get_resize(tempI2.width()*eta, tempI2.height()*eta, tempI2.depth()*eta, 1, 5);
//	//	//tempI2 = timg2;
//	//	inputVs.insert(_inputV.get_resize(twidth*t_eta, theight*t_eta, tdepth, 3, 5));
//
//	//	/*inputVs.insert(tempV);
//	//	image timg3 = tempV.get_resize(tempV.width()*eta, tempV.height()*eta, tempV.depth()*eta, 3, 5);
//	//	tempV = timg3;*/
//	//	inputY1s.insert(_inputY1.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, 5));
//	//	/*inputY1s.insert(tempY1);
//	//	image timg4 = tempY1.get_resize(tempY1.width()*eta, tempY1.height()*eta, tempY1.depth()*eta, 3, 5);
//	//	tempY1 = timg4;*/
//
//	//	inputY2s.insert(_inputY2.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, 5));
//	//	/*	inputY2s.insert(tempY2);
//	//	image timg5 = tempY2.get_resize(tempY2.width()*eta, tempY2.height()*eta, tempY2.depth()*eta, 3, 5);
//	//	tempY2 = timg5;*/
//
//	//	inputY3s.insert(_inputY3.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, 5));
//
//	//	//inputY3s.insert(tempY3);
//	//	//image timg6 = tempY3.get_resize(tempY3.width()*eta, tempY3.height()*eta, tempY3.depth()*eta, 3, 5);
//	//	//tempY3 = timg6;
//
//	//	inputY4s.insert(_inputY4.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//	//	/*inputY4s.insert(tempY4);
//	//	image timg7 = tempY4.get_resize(tempY4.width()*eta, tempY4.height()*eta, tempY4.depth()*eta, 1, 5);
//	//	tempY3 = timg6;*/
//
//	//	t_eta = t_eta*eta;
//	//}
//	//Im1s.reverse();
//	//Im2s.reverse();
//	//inputVs.reverse();
//	////inputYs.reverse();
//	//inputY1s.reverse();
//	//inputY2s.reverse();
//	//inputY3s.reverse();
//	//inputY4s.reverse();
//	////load old version completed
//
//
//	//imagelist dims;
//	int twidth = I1_blur.width();
//	int theight = I1_blur.height();
//	int tdepth = I1_blur.depth();
//	float t_eta = 1.0f;
//
//	int tdims[] = { twidth ,theight, tdepth };
//	//Im1s.insert(I1_blur.get_resize(twidth, theight, tdepth, 1, interpo));
//	//Im2s.insert(I2_blur.get_resize(twidth, theight, tdepth, 1, interpo));
//	//inputVs.insert(_inputV.get_resize(twidth, theight, tdepth, 3, interpo));
//	///*inputV2s.insert(_inputV.get_resize(twidth, theight, tdepth, 1, interpo));
//	//inputV3s.insert(_inputV.get_resize(twidth, theight, tdepth, 1, interpo));*/
//	//inputY1s.insert(_inputY1.get_resize(twidth, theight, tdepth, 3, interpo));
//	//inputY2s.insert(_inputY2.get_resize(twidth, theight, tdepth, 3, interpo));
//	//inputY3s.insert(_inputY3.get_resize(twidth, theight, tdepth, 3, interpo));
//	//inputY4s.insert(_inputY4.get_resize(twidth, theight, tdepth, 1, interpo));
//
//	//
//	/*Im1s.insert(I1_blur);
//	Im2s.insert(I2_blur);
//	inputVs.insert(_inputV);
//	inputY1s.insert(_inputY1);
//	inputY2s.insert(_inputY2);
//	inputY3s.insert(_inputY3);
//	inputY4s.insert(_inputY4);*/
////	for (int iscale = 0; iscale < scales-1; iscale++)
//	for (int iscale = 0; iscale < scales ; iscale++)
//	{
//
//	//	int tdims[] = { Im1s[iscale].width()*eta  ,Im1s[iscale].height()*eta , Im1s[iscale].depth()*eta };
//	////	int dim_in[] = { Im1s[iscale].width()  ,Im1s[iscale].height() , Im1s[iscale].depth()};
//	////	int aabb_in[] = { 0,dim_in[0]  ,0,dim_in[1]  , 0,dim_in[2] };
//	////	int dim_out[] = { Im1s[iscale].width()*eta  ,Im1s[iscale].height()*eta , Im1s[iscale].depth()*eta };
//	////	int aabb_out[] = { 0,dim_out[0]  ,0,dim_out[1]  , dim_out[2] };
//
//	////	cout << "dim_in[]  :" << dim_in[0] << " " << dim_in[1] << " " << dim_in[2] << endl;
//	////	cout << "dim_out[]  :" << dim_out[0] << " " << dim_out[1] << " " << dim_out[2] << endl;
//
//	//	Im1s.insert(Im1s[iscale].get_resize(tdims[0], tdims[1], tdims[2],1, interpo));
//	//	Im2s.insert(Im2s[iscale].get_resize(tdims[0], tdims[1], tdims[2], 1, interpo));
//	//	inputVs.insert(inputVs[iscale].get_resize(tdims[0], tdims[1], tdims[2], 3, interpo));
//	//	//inputV2s.insert(inputV2s[iscale].get_resize(tdims[0], tdims[1], tdims[2], 1, interpo));
//	//	//inputV3s.insert(inputV3s[iscale].get_resize(tdims[0], tdims[1], tdims[2], 1, interpo));
//	//	inputY1s.insert(inputY1s[iscale].get_resize(tdims[0], tdims[1], tdims[2], 3, interpo));
//	//	inputY2s.insert(inputY2s[iscale].get_resize(tdims[0], tdims[1], tdims[2], 3, interpo));
//	//	inputY3s.insert(inputY3s[iscale].get_resize(tdims[0], tdims[1], tdims[2], 3, interpo));
//	//	inputY4s.insert(inputY4s[iscale].get_resize(tdims[0], tdims[1], tdims[2], 1, interpo));
//
//		//for (std::map<std::string, image>::iterator it = F.cell_data.begin(); it != F.cell_data.end(); ++it) {
//			//image &input = it->second;
//			//if (blur > 0.001)
//			//	input.blur(blur);
//
//		//	int nc = 3;
//		//	image tmp1s(dim_out[0], dim_out[1], dim_out[2], 1);
//		//	image tmp2s(dim_out[0], dim_out[1], dim_out[2], 1);
//		//	image tmpY4s(dim_out[0], dim_out[1], dim_out[2], 1);
//		//	image tmpVs(dim_out[0], dim_out[1], dim_out[2], 3);
//		//	image tmpY1s(dim_out[0], dim_out[1], dim_out[2], 3);
//		//	image tmpY2s(dim_out[0], dim_out[1], dim_out[2], 3);
//		//	image tmpY3s(dim_out[0], dim_out[1], dim_out[2], 3);
//		//	for (int k = 0; k<dim_out[2]; k++) {
//		//		for (int j = 0; j<dim_out[1]; j++) {
//		//			for (int i = 0; i<dim_out[0]; i++) {
//		//				real3 p = advection::cell_grid_to_world(dim_out, aabb_out, real3(i, j, k));
//		//				tmp1s(i, j, k, 0) = advection::sample_image(dim_in, aabb_in, Im1s[iscale], p, 0);
//		//				tmp2s(i, j, k, 0) = advection::sample_image(dim_in, aabb_in, Im2s[iscale], p, 0);
//		//				tmpY4s(i, j, k, 0) = advection::sample_image(dim_in, aabb_in, inputY4s[iscale], p, 0);
//		//			}
//		//		}
//		//	}
//		//	Im1s.insert(tmp1s);
//		//	Im2s.insert(tmp2s);
//		//	inputY4s.insert(tmpY4s);
//
//		//	for (int c = 0; c<nc; c++) {
//		//		for (int k = 0; k<dim_out[2]; k++) {
//		//				for (int j = 0; j<dim_out[1]; j++) {
//		//					for (int i = 0; i<dim_out[0]; i++) {					
//		//					real3 p = advection::cell_grid_to_world(dim_out, aabb_out, real3(i, j, k));
//		//					tmpVs(i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputVs[iscale], p, c);
//		//					tmpY1s(i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputY1s[iscale], p, c);
//		//					tmpY2s(i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputY2s[iscale], p, c);
//		//					tmpY3s(i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputY3s[iscale], p, c);
//		//				}
//		//			}
//		//		}
//		//	}
//		//	//it->second = tmp;
//		//
//		////}	
//		//	inputVs.insert(tmpVs);
//		//	inputY1s.insert(tmpY1s);
//		//	inputY2s.insert(tmpY2s);
//		//	inputY3s.insert(tmpY3s);
//
//		//inputY1s[i].width
//
//		Im1s.insert(I1_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));
//		Im2s.insert(I2_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));
//		inputVs.insert(_inputV.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
//		inputY1s.insert(_inputY1.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
//		inputY2s.insert(_inputY2.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
//		inputY3s.insert(_inputY3.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
//		inputY4s.insert(_inputY4.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));
//
//		t_eta = t_eta*eta;
//
//
//
//
//
//
//		/*	inputY4s.insert(_inputY4.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, 5));
//		inputY4s.insert(tempY4);
//		image timg7 = tempY4.get_resize(tempY4.width()*eta, tempY4.height()*eta, tempY4.depth()*eta, 1, 5);
//		tempY3 = timg6;*/
//
//		/*Im1s.insert(tempI1);
//		image timg = tempI1.get_resize(tempI1.width()*eta, tempI1.height()*eta, tempI1.depth(), 1, 5);
//		tempI1 = timg;
//
//		Im2s.insert(tempI2);
//		image timg2 = tempI2.get_resize(tempI2.width()*eta, tempI2.height()*eta, tempI2.depth(), 1, 5);
//		tempI2 = timg2;
//
//		inputVs.insert(tempV);
//		image timg3 = tempV.get_resize(tempV.width()*eta, tempV.height()*eta, tempV.depth(), 2, 5);
//		tempV = timg3;
//
//		inputY1s.insert(tempY1);
//		image timg4 = tempY1.get_resize(tempY1.width()*eta, tempY1.height()*eta, tempY1.depth(), 2, 5);
//		tempY1 = timg4;
//
//		inputY2s.insert(tempY2);
//		image timg5 = tempY2.get_resize(tempY2.width()*eta, tempY2.height()*eta, tempY2.depth(), 2, 5);
//		tempY2 = timg5;
//
//		inputY3s.insert(tempY3);
//		image timg6 = tempY3.get_resize(tempY3.width()*eta, tempY3.height()*eta, tempY3.depth(), 1, 5);
//		tempY3 = timg6;
//		*/
//	}
//	Im1s.reverse();
//	Im2s.reverse();
//	inputVs.reverse();
//	//inputV2s.reverse();
//	//inputV3s.reverse();
//	//inputYs.reverse();
//	inputY1s.reverse();
//	inputY2s.reverse();
//	inputY3s.reverse();
//	inputY4s.reverse();
//	//cout << "Test 1" << endl;
//	////for (CImgList<>::iterator it = list.begin(); it<list.end(); ++it) (*it).mirror('x');
//	for (int iscale = 0; iscale < inputVs.size(); iscale++)
//	{
//		int levels = inputVs.size() - iscale;
//		cout << "At levels " << levels << endl;
//		int _dim[] = { Im1s[iscale].width(),Im1s[iscale].height(),Im1s[iscale].depth() };
//		//int t_dim[] = { inputY1s[iscale].width(),inputY1s[iscale].height(),inputY1s[iscale].depth() };
//		cout << "_dim[]  :" << _dim[0] << " " << _dim[1] << " " << _dim[2] << endl;
//		//cout << "_dims[]  :" << _dim[0] << " " << _dim[1] << " " << _dim[2] << endl;
//		//L1TVOpticalFlowNonlinear(_dim, Im1s(i), Im2s(i), inputVs(i), inputYs(i), _tol, _lambda, _maxIterations, _norm,
//		//									_numberOfWarps);
//
//		//real _aabb[] = { 0,_dim[0]*10000,0,_dim[1] * 10000 ,0,_dim[2] * 10000 };
//		//real _aabb[] = { _dim[0],_dim[0] * 1 ,0,_dim[1] * 1  ,0,_dim[2] * 1 };
//		//real _aabb[] = { -(_dim[0] )*0.5f,(_dim[0] )*0.5f,-(_dim[1] )*0.5f,(_dim[1])*0.5f,-(_dim[2] )*0.5f,(_dim[2] )*0.5f };
//		
//		real _aabb[] = { 0,tdims[0], 0,tdims[1], 0,tdims[2] };
//		
//			//real _aabb[] = { -tdims[0]*0.5,tdims[0]*0.5, -tdims[1]*0.5,tdims[1]*0.5, -tdims[2]*0.5,tdims[2]*0.5 };
//			cout << "tdims[]  :" << tdims[0] << " " << tdims[1] << " " << tdims[2] << endl;
//			//if ((inputVs.size() - iscale) == 1)
//			//{
//				ROI[0] *= _dim[0]/ tdims[0];
//				ROI[1] *= _dim[0]/ tdims[0];
//				ROI[2] *= _dim[1] / tdims[1];
//				ROI[3] *= _dim[1] / tdims[1];
//				ROI[4] *= _dim[2]  / tdims[2];
//				ROI[5] *= _dim[2]  / tdims[2];
//				cout << "ROI: " << ROI[0] << " " << ROI[1] << " " << ROI[2] << " " << ROI[3] << " " << ROI[4] << " " << ROI[5] << " " << endl;
//			//}
//	//	L1TVOpticalFlowNonlinear3D(_dim, Im1s[iscale], Im2s[iscale], inputV1s[iscale], inputV2s[iscale], inputV3s[iscale],
//		//	inputY1s[iscale], inputY2s[iscale],inputY3s[iscale], inputY4s[iscale], _tol, _lambda, _maxIterations, _norm,_numberOfWarps,_aabb);
//			//image test = inputVs[iscale];
//		L1TVOpticalFlowNonlinear3D(_dim, Im1s[iscale], Im2s[iscale], inputVs[iscale],
//			inputY1s[iscale], inputY2s[iscale], inputY3s[iscale], inputY4s[iscale], _tol, _lambda, _maxIterations, _norm, _numberOfWarps, _aabb,_huber);
//			//inputVs[iscale] = test;
//	
//
//		if (iscale < inputVs.size() - 1)
//		{
////		//	// use grid interpolation (4) not cubic interpolation (5)
////		//	int dim_in[] = { inputVs[iscale].width(), inputVs[iscale].height(), inputVs[iscale].depth()};
////		//	//int aabb_in[] = { inputVs[i ].width(), inputVs[i].height(), inputVs[i].depth() };
////		//	int aabb_in[] = { 0,dim_in[0]  ,0,dim_in[1]  , 0,dim_in[2] };
////		//	int dim_out[] = { inputVs[iscale + 1].width(), inputVs[iscale + 1].height(), inputVs[iscale + 1].depth() };
////		//	//int aabb_out[] = { inputVs[i + 1].width(), inputVs[i + 1].height(), inputVs[i + 1].depth() };
////		//	int aabb_out[] = { 0,dim_out[0]  ,0,dim_out[1]  , dim_out[2] };
////		
////		//	m_aabb[0] = 0 * 100;
////		//	m_aabb[1] = (_dims[0]-1) * 100;
////		//	m_aabb[2] = 0 * 100;
////		//	m_aabb[3] = (_dims[1]-1) * 100;
////		//	m_aabb[4] = 0 * 100;
////		//	m_aabb[5] = (_dims[2]-1) * 100;*/
////		//	//inputVs[i + 1]
////
////		//	image tmpY4s(dim_out[0], dim_out[1], dim_out[2], 1);
////		//	imagelist tmpVs(dim_out[0], dim_out[1], dim_out[2], 3);
////		//	imagelist tmpY1s(dim_out[0], dim_out[1], dim_out[2], 3);
////		//	imagelist tmpY2s(dim_out[0], dim_out[1], dim_out[2], 3);
////		//	imagelist tmpY3s(dim_out[0], dim_out[1], dim_out[2], 3);
////
////
////		//	// do not know if can do upsampling.
////
////		//	int nc = 3;
////		//	for (int c = 0; c<nc; c++) {
////		//		for (int k = 0; k<dim_out[2]; k++) {
////		//			for (int j = 0; j<dim_out[1]; j++) {
////		//				for (int i = 0; i<dim_out[0]; i++) {
////		//					real3 p = advection::cell_grid_to_world(dim_out, aabb_out, real3(i, j, k));
////		//					inputVs[iscale + 1](i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputVs[iscale], p, c);
////		//					inputY1s[iscale + 1](i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputY1s[iscale], p, c);
////		//					inputY2s[iscale + 1](i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputY2s[iscale], p, c);
////		//					inputY3s[iscale + 1](i, j, k, c) = advection::sample_image(dim_in, aabb_in, inputY3s[iscale], p, c);
////		//				}
////		//			}
////		//		}
////		//	}
////
////		//	for (int k = 0; k<dim_out[2]; k++) {
////		//		for (int j = 0; j<dim_out[1]; j++) {
////		//			for (int i = 0; i<dim_out[0]; i++) {
////		//				real3 p = advection::cell_grid_to_world(dim_out, aabb_out, real3(i, j, k));
////		//				inputY4s[iscale + 1](i, j, k, 0) = advection::sample_image(dim_in, aabb_in, inputY4s[iscale], p, 0);
////		//			}
////		//		}
////		//	}
////			// we can upsample but with channel each.
////		//	int _dims[] = { inputVs[iscale + 1].width(), inputVs[iscale + 1].height(), inputVs[iscale + 1].depth() };
////			
////			
////		/*	
////		int nc = 3;
////		for (size_t ic = 0; ic < nc; ic++)
////			{
////				
////				
////				inputVs[iscale + 1] = inputVs[iscale].get_resize(dim_out[0], dim_out[1], dim_out[2], ic, interpo);
////				inputY1s[iscale + 1] = inputY1s[iscale].get_resize(dim_out[0], dim_out[1], dim_out[2], ic, interpo);
////				inputY2s[iscale + 1] = inputY2s[iscale].get_resize(dim_out[0], dim_out[1], dim_out[2], ic, interpo);
////				inputY3s[iscale + 1] = inputY3s[iscale].get_resize(dim_out[0], dim_out[1], dim_out[2], ic, interpo);
////
////			}*/
////		/*	inputVs[iscale + 1] = inputVs[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
////			inputY1s[iscale + 1] = inputY1s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
////			inputY2s[iscale + 1] = inputY2s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
////			inputY3s[iscale + 1] = inputY3s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
////			inputY4s[iscale + 1] = inputY4s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 1, interpo);
////*/
////
////			//cout<< inputY1s[iscale + 1](_dims[0]*0.5, _dims[1]*, _dims[2])
//		int _dims[] = { inputVs[iscale + 1].width(), inputVs[iscale + 1].height(), inputVs[iscale + 1].depth() };
////			int _olddims[] = { inputVs[iscale].width(), inputVs[iscale].height(), inputVs[iscale].depth() };
////
////			image tmpV1(_olddims[0], _olddims[1], _olddims[2], 1, 0);
////			image tmpV2(_olddims[0], _olddims[1], _olddims[2], 1, 0);
////			image tmpV3(_olddims[0], _olddims[1], _olddims[2], 1, 0);
////			image newtmpV1(_dims[0], _dims[1], _dims[2], 1, 0);
////			image newtmpV2(_dims[0], _dims[1], _dims[2], 1, 0);
////			image newtmpV3(_dims[0], _dims[1], _dims[2], 1, 0);
////			
////			for (int k = 0; k < _olddims[2]; k++) {
////				for (int j = 0; j < _olddims[1]; j++) {
////					for (int i = 0; i < _olddims[0]; i++) {
////						//tmpVs.insert(inputVs[iscale].get_resize(_dims[0], _dims[1], _dims[2], ic, interpo));
////						tmpV1(i, j, k) = inputVs[iscale](i, j, k, 0);
////						tmpV2(i, j, k) = inputVs[iscale](i, j, k, 1);
////						tmpV3(i, j, k) = inputVs[iscale](i, j, k, 2);
////						/*if (iscale>2)
////						{
////						cout << "tmpV1(i, j, k):	" << tmpV1(i, j, k) << " " << tmpV2(i, j, k) << " " << tmpV3(i, j, k) << endl;
////						}*/
////					}
////				}
////			}
////			
////			// the center can be setted by resize fuction!!
////			newtmpV1 = tmpV1.get_resize(_dims[0], _dims[1], _dims[2],1, interpo);
////			newtmpV2 = tmpV2.get_resize(_dims[0], _dims[1], _dims[2], 1, interpo);
////			newtmpV3 = tmpV3.get_resize(_dims[0], _dims[1], _dims[2], 1, interpo);
////
////			image _V(_dims[0], _dims[1], _dims[2], 3,0) ;
////
////			for (int k = 0; k < _dims[2]; k++) {
////				for (int j = 0; j < _dims[1]; j++) {
////					for (int i = 0; i < _dims[0]; i++) {
////						//tmpVs.insert(inputVs[iscale].get_resize(_dims[0], _dims[1], _dims[2], ic, interpo));
////						_V(i, j, k, 0) = newtmpV1(i, j, k);
////						_V(i, j, k, 1) = newtmpV2(i, j, k);
////						_V(i, j, k, 2) = newtmpV3(i, j, k);
////						//cout << " _V(i, j, k) " << _V(i, j, k, 0) << " " << _V(i, j, k, 0) << " " << _V(i, j, k, 0) << endl;
////					}
////				}
////			}
////				/*image v1 = inputVs[iscale];
////				image v2 = inputVs[iscale];
////				image v3 = inputVs[iscale];
////*/
////				
////				
////			
////
////		/*	cout << "tmpVs.size()	:	"<<tmpVs.size() << endl;
////			image tmpVs2(_dims[0], _dims[1], _dims[2],3,0);
////			for (int k = 0; k < _dims[2]; k++) {
////				for (int j = 0; j < _dims[1]; j++) {
////					for (int i = 0; i < _dims[0]; i++) {
////						tmpVs2(i, j, k, 0) = tmpVs[0](i, j, k);
////						tmpVs2(i, j, k, 1) = tmpVs[1](i, j, k);
////						tmpVs2(i, j, k, 2) = tmpVs[2](i, j, k);
////
////					}
////				}
////			}*/
////
////			//inputVs[iscale + 1] = _V;
////			////inputVs[i + 1]
//			inputVs[iscale + 1] = inputVs[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
//			 
//			inputY1s[iscale + 1] = inputY1s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
//			inputY2s[iscale + 1] = inputY2s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
//			inputY3s[iscale + 1] = inputY3s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
//			inputY4s[iscale + 1] = inputY4s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 1, interpo);
//
//			//for (size_t idx = 0; idx < _dims[0]* _dims[1]* _dims[2]; idx++)
//			//{
//			/*	cout << "inputY1s: " << inputVs[iscale + 1](_dims[0] - 5, _dims[1] - 5, _dims[2] - 5, 0) << " "
//					<< inputVs[iscale + 1](_dims[0] - 5, _dims[1] - 5, _dims[2] - 5, 1) << " "
//					<< inputVs[iscale + 1](_dims[0] - 5, _dims[1] - 5, _dims[2] - 5, 2) << " " << endl;
//				cout << "inputY1s: " << inputVs[iscale + 1](_dims[0]*0.5, _dims[1] * 0.5, _dims[2] * 0.5, 0) << " "
//					<< inputVs[iscale + 1](_dims[0] * 0.5, _dims[1] * 0.5, _dims[2] * 0.5, 1) << " "
//					<< inputVs[iscale + 1](_dims[0] * 0.5, _dims[1] * 0.5, _dims[2] * 0.5, 2) << " " << endl;*/
//			//}
//				//for (int k = 0; k<_dims[2]; k++) {
//				//	for (int j = 0; j<_dims[1]; j++) {
//				//		for (int i = 0; i<_dims[0]; i++) {
//
//				//			if(inputVs[iscale + 1](i, j, k, 0)>0.0f || inputVs[iscale + 1](i, j, k, 0)<0.0f)
//				//			{ 
//				//			cout << "inputVS: " << inputVs[iscale + 1](i, j, k, 0) << " "
//				//				<< inputVs[iscale + 1](i,j,k, 1) << " "
//				//				<< inputVs[iscale + 1](i,j,k, 2) << " " << endl;
//
//				//			/*cout << "inputY2s: " << inputY2s[iscale + 1](i, j, k, 0) << " "
//				//				<< inputY2s[iscale + 1](i, j, k, 1) << " "
//				//				<< inputY2s[iscale + 1](i, j, k, 2) << " " << endl;*/
//				//			}
//				//			/*cout << "inputY1s: " << inputVs[iscale + 1](_dims[0] * 0.5, _dims[1] * 0.5, _dims[2] * 0.5, 0) << " "
//				//				<< inputVs[iscale + 1](_dims[0] * 0.5, _dims[1] * 0.5, _dims[2] * 0.5, 1) << " "
//				//				<< inputVs[iscale + 1](_dims[0] * 0.5, _dims[1] * 0.5, _dims[2] * 0.5, 2) << " " << endl;*/
//				//			//}
//				//		}
//				//	}
//				//}
//		}
//
//
//		//cout << inputVs.size()<<" "<< inputVs(i).width() << " " << inputVs(i).height() << endl;
//
//	}
//	////cout << Im1s.size() << endl;
//	////system("pause");
//
//}
//
////template< typename real, typename image >
//
//
////#if DIMENSION_SIZE ==2
//
////void doWarp(const float *image1, const float *image2, float *v1, float *v2, float *ut, float *ux, float *uy, 
////	float *uxt, float *uyt, float *uxx, float *uxy, float *uyx, float *uyy, const int *sizeMat)
//#endif
//
////
////void doWarp( float *image1,  float *image2, float *v1, float *v2, float *ut, float *ux, float *uy,  int *sizeMat)
////{
////
////	//		//doWarp(image1f, image2f, v1, v2, ut, ux, uy, sizeImage);
////
////	int nPx = (int)(sizeMat[0] * sizeMat[1]);
////    
////    float* origX = new float[nPx];
////    float* origY = new float[nPx];
////    float* origXp = new float[nPx];
////    float* origXm = new float[nPx];
////    float* origYp = new float[nPx];
////    float* origYm = new float[nPx];
////    
////	float* shiftX = new float[nPx];
////	float* shiftY = new float[nPx];
////	float* shiftXp = new float[nPx];
////	float* shiftXm = new float[nPx];
////	float* shiftYp = new float[nPx];
////	float* shiftYm = new float[nPx];
////    //create indicator grids
////	#pragma omp parallel for
////	for (int j = 0; j < sizeMat[1]; ++j)
////	{
////		for (int i = 0; i < sizeMat[0]; ++i)
////		{
////			int tmpIndex = index2DtoLinear(sizeMat, i, j);
////            
////            origX[tmpIndex] = float(i);
////            origY[tmpIndex] = float(j);
////            shiftX[tmpIndex] = float(i) + v2[tmpIndex];
////            shiftY[tmpIndex] = float(j) + v1[tmpIndex];
////            
////            shiftXp[tmpIndex] = myMin(sizeMat[0],myMax(0,shiftX[tmpIndex] + 0.5f));
////            shiftXm[tmpIndex] = myMin(sizeMat[0],myMax(0,shiftX[tmpIndex] - 0.5f));
////            shiftYp[tmpIndex] = myMin(sizeMat[1],myMax(0,shiftY[tmpIndex] + 0.5f));
////            shiftYm[tmpIndex] = myMin(sizeMat[1],myMax(0,shiftY[tmpIndex] - 0.5f));
////            
////            origXp[tmpIndex] = myMin(sizeMat[0],myMax(0,origX[tmpIndex] + 0.5f));
////            origXm[tmpIndex] = myMin(sizeMat[0],myMax(0,origX[tmpIndex] - 0.5f));
////            origYp[tmpIndex] = myMin(sizeMat[1],myMax(0,origY[tmpIndex] + 0.5f));
////            origYm[tmpIndex] = myMin(sizeMat[1],myMax(0,origY[tmpIndex] - 0.5f));
////        }
////    }
////
////    #pragma omp parallel for
////	for (int j = 0; j < sizeMat[1]; ++j)
////	{
////		for (int i = 0; i < sizeMat[0]; ++i)
////		{
////			int tmpIndex = index2DtoLinear(sizeMat, i, j);
////            
////			uy[tmpIndex] = cubicInterpolation(image2,shiftXp[tmpIndex],shiftY[tmpIndex], sizeMat) - cubicInterpolation(image2,shiftXm[tmpIndex],shiftY[tmpIndex], sizeMat);
////			ux[tmpIndex] = cubicInterpolation(image2,shiftX[tmpIndex],shiftYp[tmpIndex], sizeMat) - cubicInterpolation(image2,shiftX[tmpIndex],shiftYm[tmpIndex], sizeMat);
////			ut[tmpIndex] = cubicInterpolation(image2,shiftX[tmpIndex],shiftY[tmpIndex], sizeMat) - image1[tmpIndex] - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex];
////		}
////	}
////    
////    #pragma omp parallel for
////	for (int j = 0; j < sizeMat[1]; ++j)
////	{
////		for (int i = 0; i < sizeMat[0]; ++i)
////		{
////			int tmpIndex = index2DtoLinear(sizeMat, i, j);
////
////            float image1x = cubicInterpolation(image1,origX[tmpIndex],origYp[tmpIndex], sizeMat) - cubicInterpolation(image1,origX[tmpIndex],origYm[tmpIndex], sizeMat);
////            float image1y = cubicInterpolation(image1,origXp[tmpIndex],origY[tmpIndex], sizeMat) - cubicInterpolation(image1,origXm[tmpIndex],origY[tmpIndex], sizeMat);
////            
////		/*	uyy[tmpIndex] = cubicInterpolation(uy,origXp[tmpIndex],origY[tmpIndex], sizeMat) - cubicInterpolation(uy,origXm[tmpIndex],origY[tmpIndex], sizeMat);
////            uyx[tmpIndex] = cubicInterpolation(uy,origX[tmpIndex],origYp[tmpIndex], sizeMat) - cubicInterpolation(uy,origX[tmpIndex],origYm[tmpIndex], sizeMat);
////            uxy[tmpIndex] = cubicInterpolation(ux,origXp[tmpIndex],origY[tmpIndex], sizeMat) - cubicInterpolation(ux,origXm[tmpIndex],origY[tmpIndex], sizeMat);
////            uxx[tmpIndex] = cubicInterpolation(ux,origX[tmpIndex],origYp[tmpIndex], sizeMat) - cubicInterpolation(ux,origX[tmpIndex],origYm[tmpIndex], sizeMat);
////            */
////           // uxt[tmpIndex] = ux[tmpIndex]-image1x - uxx[tmpIndex] * v1[tmpIndex] - uxy[tmpIndex] * v2[tmpIndex];
////           // uyt[tmpIndex] = uy[tmpIndex]-image1y - uyx[tmpIndex] * v1[tmpIndex] - uyy[tmpIndex] * v2[tmpIndex];
////            
////            //check out of image:
////            //due to derivative calculation, go one pixel inside
////            if (shiftX[tmpIndex]<1.0f || shiftY[tmpIndex]<1.0f || shiftX[tmpIndex] > (sizeMat[0]-2.0f) || shiftY[tmpIndex] > (sizeMat[1]-2.0f)
////            || origX[tmpIndex]<1.0f || origY[tmpIndex]<1.0f || origX[tmpIndex] > (sizeMat[0]-2.0f) || origY[tmpIndex] > (sizeMat[1]-2.0f))
////            {
////                uy[tmpIndex] = 0.0f;
////                ux[tmpIndex] = 0.0f;
////                ut[tmpIndex] = 0.0f;
////               /* 
////                uxx[tmpIndex] = 0.0f;
////                uxy[tmpIndex] = 0.0f;
////                uyx[tmpIndex] = 0.0f;
////                uyy[tmpIndex] = 0.0f;*/
////                
////              //  uxt[tmpIndex] = 0.0f;
////              //  uyt[tmpIndex] = 0.0f;
////            }
////		}
////	}
////    
////    delete[] origX;
////    delete[] origY;
////    delete[] shiftX;
////    delete[] shiftY;
////
////    delete[] shiftXp;
////    delete[] shiftXm;
////    delete[] shiftYp;
////    delete[] shiftYm;
////
////    delete[] origXp;
////    delete[] origXm;
////    delete[] origYp;
////    delete[] origYm;
////}
////
////
////
////
//
//
////}
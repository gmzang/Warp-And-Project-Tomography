#ifndef OPTICAL_FLOW_H
#define OPTICAL_FLOW_H

#include<thread>

#define cimg_use_tiff
#define cimg_use_tif
//#include<CImg.h>
#include "CImg.h"
//#include"linear_solver.h"
namespace cl = cimg_library;
//namespace ls = linear_solver;
#include<vec3.h>

#include<options.h>

#include"cg.h"
#include"advection.h"
#include "tools.h"

#include<scope.h>
utilities_scope_defines
int interpo = 5;
#include<string_utils.h>
#define DIMENSION_SIZE 6
int ROI[] = { 5,5,5,5,5,5 };




namespace optical_flow {

	class cell_indexer {
	private:
		int m_dim[3];
	public:
		cell_indexer(const int *dim) {
			m_dim[0] = dim[0];
			m_dim[1] = dim[1];
			m_dim[2] = dim[2];
		}

		inline bool is_boundary(int i, int j, int k) const {
			return i < 0 || i >= m_dim[0] || j < 0 || j >= m_dim[1] || k < 0 || k >= m_dim[2];
		}

		inline int num_cells() const {
			return m_dim[0] * m_dim[1] * m_dim[2];
		}

		inline int operator()(int i, int j, int k) const {
			i = std::max(0, std::min(m_dim[0] - 1, i));
			j = std::max(0, std::min(m_dim[1] - 1, j));
			k = std::max(0, std::min(m_dim[2] - 1, k));
			return i + m_dim[0] * (j + k * m_dim[1]);
		}
	};

	class cell_var_indexer : public cell_indexer {
	private:
		int m_nvars;
	public:
		cell_var_indexer(const int *dim, const int nvars) : cell_indexer(dim), m_nvars(nvars) { }
		inline int operator()(int i, int j, int k, int c) const {
			return cell_indexer::operator()(i, j, k)*m_nvars + c;
		}
		inline int num_vars() {
			return m_nvars * num_cells();
		}
	};


	template< typename real, typename image >
	image blur_and_downsample(const image &I, const int nx, const int ny, const int nz, const real sigma) {
		return I.get_blur(sigma).get_resize(nx, ny, nz, I.spectrum(), 5);
	}

	template< typename real, typename image >
	image blur_and_downsample(const image &I, const real eta, const real sigma) {
		return blur_and_downsample(I, I.width()*eta, I.height()*eta, I.depth(), sigma);
	}
	template< typename real, typename image >
	image Y_blur_and_downsample(const image &I, const real eta, const real sigma) {
		return blur_and_downsample(I, I.width()*eta, I.height()*eta, I.depth(), sigma);
	}

	template< typename real, typename image >
	void L1TVOpticalFlowNonlinear3D(int *dim, const image & image1, const image & image2, image & _inputV,
		image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4,
		float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps, real *_aabb, float huber, int *m_roi)
		//void L1TVOpticalFlowNonlinear3D(int *dim, const image & image1, const image & image2, image & _inputV1, image & _inputV2, image & _inputV3, 
		//	image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4,
		//	float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps, real *_aabb)
	{

		//const int maxIterations = (int)mxGetScalar(prhs[4]);

		//const size_t *sizeImage = mxGetDimensions(prhs[0]);

		//float *u1 = image1;
		//float *u2 = image2;
		//image u1 = image1;
		//image u2 = image2;

		cout << "AABB: " << _aabb[0] << " " << _aabb[1] << " " << _aabb[2] << " " << _aabb[3] << " " << _aabb[4] << " " << _aabb[5] << " " << endl;
		const float tol = _tol;
		const float lambda = _lambda;

		//mexPrintf("tol: %f, lambda: %f; ", tol, lambda);

		const int maxIterations = _maxIterations;

		//const size_t *sizeImage = mxGetDimensions(prhs[0]);
		//int *sizeImage = new int[2];
		//for (int i = 0; i < _dims; i++)
		//{
		//	sizeImage[i] = imsizeX;
		//	sizeImage[i] = imsizeX;
		//}
		//sizeImage[0] = dim[0];
		//sizeImage[1] = dim[1];
		cout << "image1 " << image1.width() << " " << image1.height() << " " << image1.depth() << endl;
		cout << "image2 " << image2.width() << " " << image2.height() << " " << image2.depth() << endl;
		cout << "sizeImage " << dim[0] << " " << dim[1] << " " << dim[2] << endl;
		//cout << "Input " << _inputV(30,30,0,0) << " " << _inputV(80, 80, 0, 0) << endl;
		cout << " in the L1TVOpticalFlowNonlinear " << endl;

		/*float *inputV = 0;
		float *inputY = 0;*/
		//image  inputV = _inputV;
		//image  inputY = _inputY;

		int typeNorm = _norm;


		float stepsize[4] = { 1.0f, 1.0f, 1.0f ,1.0f };
		//float stepsize[3] = { _stepsize[0], _stepsize[1], _stepsize[2] };

		float stepsizeD[4] = { 1.0f / stepsize[0],
			1.0f / stepsize[1],
			1.0f / stepsize[2],
			1.0f / stepsize[3] };

		int numberOfWarps = _numberOfWarps;
		//if (nrhs > 10)
		//{
		//	numberOfWarps = (int)mxGetScalar(prhs[10]);
		//}

		//float huberEpsilon = 0.01f;
		float huberEpsilon = huber;

		//float huberEpsilon = 0.1f;
		//if (nrhs > 11)
		//{
		//	huberEpsilon = (float)mxGetScalar(prhs[11]);
		//}

		int gradientConstancy = 0;
		//if (nrhs > 12)
		//{
		//	gradientConstancy = (int)mxGetScalar(prhs[12]);
		//}
		//
		const int nPx = (int)(dim[0] * dim[1] * dim[2]);
		cout << "nPX" << nPx << endl;
		


		float* v1 = new float[nPx];
		float* v2 = new float[nPx];
		float* v3 = new float[nPx];

		float* v1Old = new float[nPx];
		float* v2Old = new float[nPx];
		float* v3Old = new float[nPx];

		float* image1f = new float[nPx];
		float* image2f = new float[nPx];

		float* ux = new float[nPx];
		float* uy = new float[nPx];
		float* uz = new float[nPx];
		float* ut = new float[nPx];
		//	image uXYZ(dim[0], dim[1], dim[2], 2, 0);
		/*float* uxx = new float[nPx];
		float* uxy = new float[nPx];
		float* uyx = new float[nPx];
		float* uyy = new float[nPx];*/

		//	float* uxt = new float[nPx];
		//float* uyt = new float[nPx];

		float* y11 = new float[nPx];
		float* y12 = new float[nPx];
		float* y13 = new float[nPx];

		float* y21 = new float[nPx];
		float* y22 = new float[nPx];
		float* y23 = new float[nPx];

		float* y31 = new float[nPx];
		float* y32 = new float[nPx];
		float* y33 = new float[nPx];
		float* y4 = new float[nPx];
		/*image Y1(dim[0], dim[1], dim[2], 2, 0);
		image Y2(dim[0], dim[1], dim[2], 2, 0);
		image Y3(dim[0], dim[1], dim[2], 1, 0);*/
		//float* y4 = new float[nPx];
		//float* y5 = new float[nPx];

		float* y11Old = new float[nPx];
		float* y12Old = new float[nPx];
		float* y13Old = new float[nPx];

		float* y21Old = new float[nPx];
		float* y22Old = new float[nPx];
		float* y23Old = new float[nPx];

		float* y31Old = new float[nPx];
		float* y32Old = new float[nPx];
		float* y33Old = new float[nPx];
		//float* y4Old = new float[nPx];


		/*image Y1old(dim[0], dim[1], dim[2], 2, 0);
		image Y2old(dim[0], dim[1], dim[2], 2, 0);
		image Y3old(dim[0], dim[1], dim[2], 1, 0);*/
		float* y4Old = new float[nPx];
		float* y5Old = new float[nPx];

		float* Kty1 = new float[nPx];
		float* Kty2 = new float[nPx];
		float* Kty3 = new float[nPx];

		float* Kty1Old = new float[nPx];
		float* Kty2Old = new float[nPx];
		float* Kty3Old = new float[nPx];

		float* Kx11 = new float[nPx];
		float* Kx12 = new float[nPx];
		float* Kx13 = new float[nPx];

		float* Kx21 = new float[nPx];
		float* Kx22 = new float[nPx];
		float* Kx23 = new float[nPx];

		float* Kx31 = new float[nPx];
		float* Kx32 = new float[nPx];
		float* Kx33 = new float[nPx];
		//float* Kx4 = new float[nPx];
		image KX1(dim[0], dim[1], dim[2], 2, 0);
		image KX2(dim[0], dim[1], dim[2], 2, 0);
		image KX3(dim[0], dim[1], dim[2], 1, 0);
		float* Kx4 = new float[nPx];
		float* Kx5 = new float[nPx];

		float* Kx11Old = new float[nPx];
		float* Kx12Old = new float[nPx];
		float* Kx13Old = new float[nPx];

		float* Kx21Old = new float[nPx];
		float* Kx22Old = new float[nPx];
		float* Kx23Old = new float[nPx];

		float* Kx31Old = new float[nPx];
		float* Kx32Old = new float[nPx];
		float* Kx33Old = new float[nPx];
		//float* Kx4Old = new float[nPx];

		//image KX1old(dim[0], dim[1], dim[2], 2, 0);
		//image KX2old(dim[0], dim[1], dim[2], 2, 0);
		//image KX3old(dim[0], dim[1], dim[2], 1, 0);
		float* Kx4Old = new float[nPx];
		float* Kx5Old = new float[nPx];

		//float sigma1 = myMin(stepsize[1] / 2.0f, stepsize[2] / 2.0f);

		float sigma1 = myMin(stepsize[3] / 3.0f, myMin(stepsize[1] / 3.0f, stepsize[2] / 3.0f));
		float* sigma2 = new float[nPx];

		//for gradient constancy
		//float* sigma3 = new float[nPx];
		//float* sigma4 = new float[nPx];

		float* tau1 = new float[nPx];
		float* tau2 = new float[nPx];
		float* tau3 = new float[nPx];

		int * tableI = new int[nPx];
		int * tableJ = new int[nPx];
		int * tableK = new int[nPx];

		//Huber Factor

		const float huberFactor = 1.0f / (1.0f + sigma1* huberEpsilon / lambda);
		//const float huberFactor = 1.0;
		//residuals
		float p = 0.0f;
		float d = 0.0f;

		#pragma omp parallel for
		for (int k = 0; k < dim[2]; ++k)
		{
			for (int j = 0; j < dim[1]; ++j)
			{
				for (int i = 0; i < dim[0]; ++i)
				{
					//int tmpIndex = index2DtoLinear(dim, i, j);
					int tmpIndex = index3DtoLinear(dim, i, j, k);

					tableI[tmpIndex] = i;
					tableJ[tmpIndex] = j;
					tableK[tmpIndex] = k;

				
					/*v1[tmpIndex] = (float)inputV[tmpIndex];
					v2[tmpIndex] = (float)inputV[nPx + tmpIndex];*/

					/*		v1[tmpIndex] = (float)inputV(tmpIndex);
					v2[tmpIndex] = (float)inputV(nPx + tmpIndex);*/
					v1[tmpIndex] = _inputV(i, j, k, 0);
					v2[tmpIndex] = _inputV(i, j, k, 1);
					v3[tmpIndex] = _inputV(i, j, k, 2);
					//_inputY1(i, j, 0, 0);
					//v1[tmpIndex] = 0.0f;
					//v2[tmpIndex] = 0.0f;


					//cout << "v1[tmpIndex]: " << v1[tmpIndex] << " " << v2[tmpIndex] << " " << v3[tmpIndex] << endl;
					Kty1[tmpIndex] = 0.0f;
					Kty2[tmpIndex] = 0.0f;
					Kty3[tmpIndex] = 0.0f;

					//if (nrhs > 7)
					//{
					//	y11[tmpIndex] = (float)inputY[0 * nPx + tmpIndex];
					//	y12[tmpIndex] = (float)inputY[1 * nPx + tmpIndex];
					//	y21[tmpIndex] = (float)inputY[2 * nPx + tmpIndex];
					//	y22[tmpIndex] = (float)inputY[3 * nPx + tmpIndex];
					//             y3[tmpIndex] = (float)inputY[4 * nPx + tmpIndex];
					//             
					//             //for gradient constancy
					//             y4[tmpIndex] = (float)inputY[5 * nPx + tmpIndex];
					//             y5[tmpIndex] = (float)inputY[6 * nPx + tmpIndex];
					//}
					//else
					//{
					//	y11[tmpIndex] = 0.0f;
					//	y12[tmpIndex] = 0.0f;
					//	y21[tmpIndex] = 0.0f;
					//	y22[tmpIndex] = 0.0f;
					//             y3[tmpIndex] = 0.0f;
					//             
					//             //for gradient constancy
					//             y4[tmpIndex] = 0.0f;
					//             y5[tmpIndex] = 0.0f;
					//}

					/*	y11[tmpIndex] = (float)_inputY(0 * nPx + tmpIndex);
					y12[tmpIndex] = (float)_inputY(1 * nPx + tmpIndex);
					y21[tmpIndex] = (float)_inputY(2 * nPx + tmpIndex);
					y22[tmpIndex] = (float)_inputY(3 * nPx + tmpIndex);
					y3[tmpIndex] = (float)inputY(4 * nPx + tmpIndex);
					*/
					y11[tmpIndex] = _inputY1(i, j, k, 0);
					y12[tmpIndex] = _inputY1(i, j, k, 1);
					y13[tmpIndex] = _inputY1(i, j, k, 2);

					y21[tmpIndex] = _inputY2(i, j, k, 0);
					y22[tmpIndex] = _inputY2(i, j, k, 1);
					y23[tmpIndex] = _inputY2(i, j, k, 2);

					y31[tmpIndex] = _inputY3(i, j, k, 0);
					y32[tmpIndex] = _inputY3(i, j, k, 1);
					y33[tmpIndex] = _inputY3(i, j, k, 2);

					y4[tmpIndex] = _inputY4(i, j, k);

					//for gradient constancy
					//	y4[tmpIndex] = (float)inputY(5 * nPx + tmpIndex);
					//	y5[tmpIndex] = (float)inputY(6 * nPx + tmpIndex);



					Kx11[tmpIndex] = 0.0f;
					Kx12[tmpIndex] = 0.0f;
					Kx13[tmpIndex] = 0.0f;

					Kx21[tmpIndex] = 0.0f;
					Kx22[tmpIndex] = 0.0f;
					Kx23[tmpIndex] = 0.0f;

					Kx31[tmpIndex] = 0.0f;
					Kx32[tmpIndex] = 0.0f;
					Kx33[tmpIndex] = 0.0f;

					Kx4[tmpIndex] = 0.0f;

					////for gradient constancy
					//Kx4[tmpIndex] = 0.0f;
					//Kx5[tmpIndex] = 0.0f;
				}
			}
		}
		cout << "Finish the initialization, open clock " << endl;
		clock_t begin = clock();
		//do k warpings
		for (int idxwarp = 0; idxwarp < numberOfWarps; ++idxwarp)
		{


			

		
			real dt = 1;
			//image _warpI2 = advection::advect(dim, _inputV, image2, -dt);
			image _warpI2 = advection::advect(dim, _aabb, _inputV, image2, -dt);
			//image _warpI2 = advection::advect(dim, _aabb, _inputV1, _inputV2, _inputV3, image2, -dt);
			imagelist warpI2_xyz = _warpI2.get_gradient("xyz", 0);
			////	image im_ux = warpI2_xy[0];
			////	image im_uy = warpI2_xy[1];

		

		#pragma omp parallel for
			for (int k = m_roi[4]; k < dim[2] - m_roi[5]; ++k)
			{
				for (int j = m_roi[2]; j < dim[1] - m_roi[3]; ++j)
				{
					for (int i = m_roi[0]; i < dim[0] - m_roi[1]; ++i)
					{
						int tmpIndex = index3DtoLinear(dim, i, j, k);
					
						ux[tmpIndex] = warpI2_xyz[0](i, j, k);
						uy[tmpIndex] = warpI2_xyz[1](i, j, k);
						uz[tmpIndex] = warpI2_xyz[2](i, j, k);
						//cout << ux[tmpIndex] << " " << uy[tmpIndex] << " " << uz[tmpIndex] << endl;
						//ut[tmpIndex] = _warpI2(i, j) - image1(i,j) - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex];
						ut[tmpIndex] = _warpI2(i, j, k) - image1(i, j, k) - ux[tmpIndex] * v1[tmpIndex] - uy[tmpIndex] * v2[tmpIndex] - uz[tmpIndex] * v3[tmpIndex];
						//ut[tmpIndex] = _warpI2(i, j) - image1(i, j) - uXYZ(i, j, 0, 0) * v1[tmpIndex] - uXYZ(i, j, 0, 1) * v2[tmpIndex];

					}
				}
			}

			//cout << " in the L1TVOpticalFlowNonlinear  2 " << endl;





			#pragma omp parallel for
			for (int i = 0; i < nPx; ++i)
			{
				//adjust step sizes
				// GM Testing 1
				tau1[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(ux[i]);
				tau2[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(uy[i]);
				tau3[i] = 4.0f / std::min(std::min(stepsize[1], stepsize[2]), stepsize[3]) + myAbs(uz[i]);
			
				//sigma2[i] = std::abs(ux[i]) + std::abs(uXYZ(i + nPx));
				sigma2[i] = std::abs(ux[i]) + std::abs(uy[i]) + std::abs(uz[i]);
				//sigma2[i] = std::abs(uXYZ(i)) + std::abs(uXYZ(i+ nPx));


				tau1[i] = 1.0f / tau1[i];
				tau2[i] = 1.0f / tau2[i];
				tau3[i] = 1.0f / tau3[i];
				sigma2[i] = 1.0f / sigma2[i];

			}


			int iterations = 0;
			float err = 1.0f;

			while (err > tol && iterations <= maxIterations)
			{
				++iterations;

				if (iterations % 50 == 0)
				{
					p = 0.0f;
					d = 0.0f;
				}

		
#pragma omp parallel for reduction(+:p)
				for (int k = m_roi[4]; k < dim[2] - m_roi[5]; ++k)
				{
					for (int j = m_roi[2]; j < dim[1] - m_roi[3]; ++j)
					{
						for (int i = m_roi[0]; i < dim[0] - m_roi[1]; ++i)
						{
							int tmpIndex = index3DtoLinear(dim, i, j, k);

						
							if (iterations % 50 == 0)
							{
								v1Old[tmpIndex] = v1[tmpIndex];
								v2Old[tmpIndex] = v2[tmpIndex];
								v3Old[tmpIndex] = v3[tmpIndex];

								Kty1Old[tmpIndex] = Kty1[tmpIndex];
								Kty2Old[tmpIndex] = Kty2[tmpIndex];
								Kty3Old[tmpIndex] = Kty3[tmpIndex];
							}


							//transpose equals -div  // m is mb  not choose 1
							Kty1[tmpIndex] = -stepsizeD[1] * dxm3(y11, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[2] * dym3(y12, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[3] * dzm3(y13, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
								ux[tmpIndex] * y4[tmpIndex];

							Kty2[tmpIndex] = -stepsizeD[1] * dxm3(y21, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[2] * dym3(y22, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[3] * dzm3(y23, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
								uy[tmpIndex] * y4[tmpIndex];

							Kty3[tmpIndex] = -stepsizeD[1] * dxm3(y31, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[2] * dym3(y32, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) -
								stepsizeD[3] * dzm3(y33, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]) +
								uz[tmpIndex] * y4[tmpIndex];

							v1[tmpIndex] -= tau1[tmpIndex] * Kty1[tmpIndex];
							v2[tmpIndex] -= tau2[tmpIndex] * Kty2[tmpIndex];
							v3[tmpIndex] -= tau3[tmpIndex] * Kty3[tmpIndex];
							_inputV(i, j, k, 0) = v1[tmpIndex];// -= tau1[i] * Kty1[i];
							_inputV(i, j, k, 1) = v2[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];
							_inputV(i, j, k, 2) = v3[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];
															   //if(v3[tmpIndex]>0.0f || v3[tmpIndex]<0.0f)
															   //	cout << v1[tmpIndex] << " " << v2[tmpIndex] << " " << v3[tmpIndex] << " " << endl;
															   //_inputV(tmpIndex) = v1[tmpIndex];// -= tau1[i] * Kty1[i];
															   //_inputV(tmpIndex + nPx) = v2[tmpIndex];//v2[i] -= tau2[i] * Kty2[i];


															   //_inputV(tmpIndex) -= tau1[tmpIndex] * Kty1[tmpIndex];
															   //_inputV(tmpIndex+ nPx) -= tau2[tmpIndex] * Kty2[tmpIndex];


							if (iterations % 50 == 0)
							{
								//residuals
								p += std::abs((v1Old[tmpIndex] - v1[tmpIndex]) / tau1[tmpIndex] - Kty1Old[tmpIndex] + Kty1[tmpIndex])
									+ std::abs((v2Old[tmpIndex] - v2[tmpIndex]) / tau2[tmpIndex] - Kty2Old[tmpIndex] + Kty2[tmpIndex])
									+ std::abs((v3Old[tmpIndex] - v3[tmpIndex]) / tau3[tmpIndex] - Kty3Old[tmpIndex] + Kty3[tmpIndex]);

							}
						}
					}
				}
				//}		//cout << " in the L1TVOpticalFlowNonlinear  4  //dual step" << endl;


				//dual step

				//imagelist V1grad=_inputV.get_shared_channel(0).get_gradient("xy", 1);
				//imagelist V2grad = _inputV.get_shared_channel(1).get_gradient("xy", 1);
				/*imagelist V1grad = _inputV.get_shared_channel(0).get_gradient("xy", 1);
				imagelist V2grad = _inputV.get_shared_channel(1).get_gradient("xy", 1);*/

				//imagelist Vgrad = _inputV.get_gradient("xy", 1);
				//image im1 = _inputV

#pragma omp parallel for reduction(+:d)
				for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
				{
					float y11Tilde, y12Tilde, y13Tilde, y21Tilde, y22Tilde, y23Tilde, y31Tilde, y32Tilde, y33Tilde;
					/* float y11Old, y12Old, y13Old, y21Old, y22Old, y23Old, y31Old, y32Old, y33Old, y4Old;
					float Kx11Old, Kx12Old, Kx13Old, Kx21Old, Kx22Old, Kx23Old, Kx31Old, Kx32Old, Kx33Old, Kx4Old;*/

		

					if (iterations % 50 == 0)
					{


						/*Y1old(tmpIndex) = _inputY1(tmpIndex);
						Y1old(tmpIndex+ nPx) = _inputY1(tmpIndex+ nPx);

						Y2old(tmpIndex) = _inputY2(tmpIndex);
						Y2old(tmpIndex + nPx) = _inputY2(tmpIndex + nPx);
						Y3old(tmpIndex) = _inputY3(tmpIndex);*/
						y11Old[tmpIndex] = y11[tmpIndex];
						y12Old[tmpIndex] = y12[tmpIndex];
						y13Old[tmpIndex] = y13[tmpIndex];
						y21Old[tmpIndex] = y21[tmpIndex];
						y22Old[tmpIndex] = y22[tmpIndex];
						y23Old[tmpIndex] = y23[tmpIndex];
						y31Old[tmpIndex] = y31[tmpIndex];
						y32Old[tmpIndex] = y32[tmpIndex];
						y33Old[tmpIndex] = y33[tmpIndex];

						y4Old[tmpIndex] = y4[tmpIndex];

				/*		y11Old = y11[tmpIndex];
						y12Old = y12[tmpIndex];
						y13Old = y13[tmpIndex];
						y21Old = y21[tmpIndex];
						y22Old = y22[tmpIndex];
						y23Old = y23[tmpIndex];
						y31Old = y31[tmpIndex];
						y32Old = y32[tmpIndex];
						y33Old = y33[tmpIndex];

						y4Old = y4[tmpIndex];*/
					}

					/*KX1old(tmpIndex) = KX1(tmpIndex);
					KX1old(tmpIndex+ nPx) = KX1(tmpIndex+ nPx);
					KX2old(tmpIndex) = KX2(tmpIndex);
					KX2old(tmpIndex + nPx) = KX2(tmpIndex + nPx);
					KX3old(tmpIndex) = KX3(tmpIndex);*/
			/*		Kx11Old = Kx11[tmpIndex];
					Kx12Old = Kx12[tmpIndex];
					Kx13Old = Kx13[tmpIndex];
					Kx21Old = Kx21[tmpIndex];
					Kx22Old = Kx22[tmpIndex];
					Kx23Old = Kx23[tmpIndex];
					Kx31Old = Kx31[tmpIndex];
					Kx32Old = Kx32[tmpIndex];
					Kx33Old = Kx33[tmpIndex];

					Kx4Old = Kx4[tmpIndex];*/

					//KX1old(tmpIndex) = KX1(tmpIndex);
					//KX1old(tmpIndex + nPx) = KX1(tmpIndex + nPx);
					//KX2old(tmpIndex) = KX2(tmpIndex);
					//KX2old(tmpIndex + nPx) = KX2(tmpIndex + nPx);
					//KX3old(tmpIndex) = KX3(tmpIndex);
					Kx11Old[tmpIndex] = Kx11[tmpIndex];
					Kx12Old[tmpIndex] = Kx12[tmpIndex];
					Kx13Old[tmpIndex] = Kx13[tmpIndex];
					Kx21Old[tmpIndex] = Kx21[tmpIndex];
					Kx22Old[tmpIndex] = Kx22[tmpIndex];
					Kx23Old[tmpIndex] = Kx23[tmpIndex];
					Kx31Old[tmpIndex] = Kx31[tmpIndex];
					Kx32Old[tmpIndex] = Kx32[tmpIndex];
					Kx33Old[tmpIndex] = Kx33[tmpIndex];

					Kx4Old[tmpIndex] = Kx4[tmpIndex];

					/*
					KX1(tmpIndex) = stepsizeD[1] * V1grad[0](tmpIndex);
					KX1(tmpIndex + nPx) = stepsizeD[2] * V1grad[1](tmpIndex);
					KX2(tmpIndex) = stepsizeD[1] * V2grad[0](tmpIndex);
					KX2(tmpIndex + nPx) = stepsizeD[2] * V2grad[1](tmpIndex);*/
					/*Kx11[tmpIndex] = stepsizeD[1] * V1grad[0](tmpIndex);
					Kx12[tmpIndex] = stepsizeD[2] * V1grad[1](tmpIndex);
					Kx21[tmpIndex] = stepsizeD[1] * V2grad[0](tmpIndex);
					Kx22[tmpIndex] = stepsizeD[2] * V2grad[1](tmpIndex);*/
					/*	cout << Kx11[tmpIndex] << " " << Kx12[tmpIndex] << " " << Kx13[tmpIndex] << " " << Kx21[tmpIndex] << " " << Kx22[tmpIndex] << " " <<
					Kx23[tmpIndex] << " " << endl;*/
					Kx11[tmpIndex] = stepsizeD[1] * dxp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx12[tmpIndex] = stepsizeD[2] * dyp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx13[tmpIndex] = stepsizeD[3] * dzp3(v1, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);

					Kx21[tmpIndex] = stepsizeD[1] * dxp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx22[tmpIndex] = stepsizeD[2] * dyp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx23[tmpIndex] = stepsizeD[3] * dzp3(v2, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);

					Kx31[tmpIndex] = stepsizeD[1] * dxp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx32[tmpIndex] = stepsizeD[2] * dyp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);
					Kx33[tmpIndex] = stepsizeD[3] * dzp3(v3, dim, tableI[tmpIndex], tableJ[tmpIndex], tableK[tmpIndex]);

					Kx4[tmpIndex] = ux[tmpIndex] * v1[tmpIndex] + uy[tmpIndex] * v2[tmpIndex] + uz[tmpIndex] * v3[tmpIndex];  // warpped_u2=(ux, uy)
																															  /*cout << Kx11[tmpIndex] << " " << Kx12[tmpIndex] << " " << Kx13[tmpIndex] << " " << Kx21[tmpIndex] << " " << Kx22[tmpIndex] << " " <<
																															  Kx23[tmpIndex] << " " << endl;*/
																															  //KX3(tmpIndex) = uXYZ(tmpIndex) * v1[tmpIndex] + uXYZ(tmpIndex + nPx) * v2[tmpIndex];  // warpped_u2=(ux, uy)
																															  //Kx3[tmpIndex] = uXYZ(tmpIndex) * v1[tmpIndex] + uXYZ(tmpIndex+ nPx) * v2[tmpIndex];  // warpped_u2=(ux, uy)

																															  /*Kx11[tmpIndex] = stepsizeD[1] * dxp(_inputV.get_channel(0), dim, tableI[tmpIndex], tableJ[tmpIndex]);
																															  Kx12[tmpIndex] = stepsizeD[2] * dyp(_inputV.get_channel(0), dim, tableI[tmpIndex], tableJ[tmpIndex]);
																															  Kx21[tmpIndex] = stepsizeD[1] * dxp(_inputV.get_channel(1), dim, tableI[tmpIndex], tableJ[tmpIndex]);
																															  Kx22[tmpIndex] = stepsizeD[2] * dyp(_inputV.get_channel(1), dim, tableI[tmpIndex], tableJ[tmpIndex]);*/

																															  //Kx3[tmpIndex] = ux[tmpIndex] * _inputV(tmpIndex) + uy[tmpIndex] * _inputV(tmpIndex+ nPx);  // warpped_u2=(ux, uy)

																															  /*if (gradientConstancy>0)
																															  {
																															  Kx4Old[tmpIndex] = Kx4[tmpIndex];
																															  Kx5Old[tmpIndex] = Kx5[tmpIndex];

																															  Kx4[tmpIndex] = uxx[tmpIndex] * v1[tmpIndex] + uxy[tmpIndex] * v2[tmpIndex];
																															  Kx5[tmpIndex] = uyx[tmpIndex] * v1[tmpIndex] + uyy[tmpIndex] * v2[tmpIndex];

																															  y4Old[tmpIndex] = y4[tmpIndex];
																															  y5Old[tmpIndex] = y5[tmpIndex];

																															  y4[tmpIndex] = std::max(-1.0f,std::min(1.0f, y4[tmpIndex] + sigma3[tmpIndex]*(Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old[tmpIndex] + uxt[tmpIndex])));
																															  y5[tmpIndex] = std::max(-1.0f,std::min(1.0f, y5[tmpIndex] + sigma4[tmpIndex]*(Kx5[tmpIndex] + Kx5[tmpIndex] - Kx5Old[tmpIndex] + uyt[tmpIndex])));
																															  }*/

					if (typeNorm == 4) // Huber
					{
						/*y11Tilde = (_inputY1(tmpIndex) + sigma1*(KX1(tmpIndex) + KX1(tmpIndex) - KX1old(tmpIndex))) * huberFactor;
						y12Tilde = (_inputY1(tmpIndex + nPx) + sigma1*(KX1(tmpIndex + nPx) + KX1(tmpIndex + nPx) - KX1old(tmpIndex + nPx))) * huberFactor;
						y21Tilde = (_inputY2(tmpIndex)+ sigma1*(KX2(tmpIndex) + KX2(tmpIndex) - KX2old(tmpIndex))) * huberFactor;
						y22Tilde = (_inputY2(tmpIndex + nPx) + sigma1*(KX2(tmpIndex + nPx) + KX2(tmpIndex + nPx) - KX2old(tmpIndex + nPx))) * huberFactor;*/
					/*	y11Tilde = (y11[tmpIndex] + sigma1 * (Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old)) * huberFactor;
						y12Tilde = (y12[tmpIndex] + sigma1 * (Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old)) * huberFactor;
						y13Tilde = (y13[tmpIndex] + sigma1 * (Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old)) * huberFactor;

						y21Tilde = (y21[tmpIndex] + sigma1 * (Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old)) * huberFactor;
						y22Tilde = (y22[tmpIndex] + sigma1 * (Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old)) * huberFactor;
						y23Tilde = (y23[tmpIndex] + sigma1 * (Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old)) * huberFactor;

						y31Tilde = (y31[tmpIndex] + sigma1 * (Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old)) * huberFactor;
						y32Tilde = (y32[tmpIndex] + sigma1 * (Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old)) * huberFactor;
						y33Tilde = (y33[tmpIndex] + sigma1 * (Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old)) * huberFactor;*/

						y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex])) * huberFactor;
						y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex])) * huberFactor;
						y13Tilde = (y13[tmpIndex] + sigma1*(Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old[tmpIndex])) * huberFactor;

						y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex])) * huberFactor;
						y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex])) * huberFactor;
						y23Tilde = (y23[tmpIndex] + sigma1*(Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old[tmpIndex])) * huberFactor;

						y31Tilde = (y31[tmpIndex] + sigma1*(Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old[tmpIndex])) * huberFactor;
						y32Tilde = (y32[tmpIndex] + sigma1*(Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old[tmpIndex])) * huberFactor;
						y33Tilde = (y33[tmpIndex] + sigma1*(Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old[tmpIndex])) * huberFactor;
					}
					else
					{
						/*y11Tilde = _inputY1(tmpIndex) + sigma1*(KX1(tmpIndex) + KX1(tmpIndex) - KX1old(tmpIndex));
						y12Tilde = _inputY1(tmpIndex + nPx) + sigma1*(KX1(tmpIndex + nPx) + KX1(tmpIndex + nPx) - KX1old(tmpIndex + nPx));
						y21Tilde = _inputY2(tmpIndex) + sigma1*(KX2(tmpIndex) + KX2(tmpIndex) - KX2old(tmpIndex));
						y22Tilde = _inputY2(tmpIndex + nPx) + sigma1*(KX2(tmpIndex + nPx) + KX2(tmpIndex + nPx) - KX2old(tmpIndex + nPx)) ;*/
						// ATV or ITV?

				/*		y11Tilde = (y11[tmpIndex] + sigma1 * (Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old));
						y12Tilde = (y12[tmpIndex] + sigma1 * (Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old));
						y13Tilde = (y13[tmpIndex] + sigma1 * (Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old));

						y21Tilde = (y21[tmpIndex] + sigma1 * (Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old));
						y22Tilde = (y22[tmpIndex] + sigma1 * (Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old));
						y23Tilde = (y23[tmpIndex] + sigma1 * (Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old));

						y31Tilde = (y31[tmpIndex] + sigma1 * (Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old));
						y32Tilde = (y32[tmpIndex] + sigma1 * (Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old));
						y33Tilde = (y33[tmpIndex] + sigma1 * (Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old));*/

						y11Tilde = (y11[tmpIndex] + sigma1*(Kx11[tmpIndex] + Kx11[tmpIndex] - Kx11Old[tmpIndex]));
						y12Tilde = (y12[tmpIndex] + sigma1*(Kx12[tmpIndex] + Kx12[tmpIndex] - Kx12Old[tmpIndex]));
						y13Tilde = (y13[tmpIndex] + sigma1*(Kx13[tmpIndex] + Kx13[tmpIndex] - Kx13Old[tmpIndex]));

						y21Tilde = (y21[tmpIndex] + sigma1*(Kx21[tmpIndex] + Kx21[tmpIndex] - Kx21Old[tmpIndex]));
						y22Tilde = (y22[tmpIndex] + sigma1*(Kx22[tmpIndex] + Kx22[tmpIndex] - Kx22Old[tmpIndex]));
						y23Tilde = (y23[tmpIndex] + sigma1*(Kx23[tmpIndex] + Kx23[tmpIndex] - Kx23Old[tmpIndex]));

						y31Tilde = (y31[tmpIndex] + sigma1*(Kx31[tmpIndex] + Kx31[tmpIndex] - Kx31Old[tmpIndex]));
						y32Tilde = (y32[tmpIndex] + sigma1*(Kx32[tmpIndex] + Kx32[tmpIndex] - Kx32Old[tmpIndex]));
						y33Tilde = (y33[tmpIndex] + sigma1*(Kx33[tmpIndex] + Kx33[tmpIndex] - Kx33Old[tmpIndex]));

					}

					float divisor1 = std::max(1.0f, std::sqrt(y11Tilde*y11Tilde + y12Tilde * y12Tilde + y13Tilde * y13Tilde) / lambda);
					float divisor2 = std::max(1.0f, std::sqrt(y21Tilde*y21Tilde + y22Tilde * y22Tilde + y23Tilde * y23Tilde) / lambda);
					float divisor3 = std::max(1.0f, std::sqrt(y31Tilde*y31Tilde + y32Tilde * y32Tilde + y33Tilde * y33Tilde) / lambda);


					

					y11[tmpIndex] = y11Tilde * lambda / (lambda + 2 * sigma1);
					y12[tmpIndex] = y12Tilde * lambda / (lambda + 2 * sigma1);
					y13[tmpIndex] = y13Tilde * lambda / (lambda + 2 * sigma1);

					y21[tmpIndex] = y21Tilde * lambda / (lambda + 2 * sigma1);
					y22[tmpIndex] = y22Tilde * lambda / (lambda + 2 * sigma1);
					y23[tmpIndex] = y23Tilde * lambda / (lambda + 2 * sigma1);

					y31[tmpIndex] = y31Tilde * lambda / (lambda + 2 * sigma1);
					y32[tmpIndex] = y32Tilde * lambda / (lambda + 2 * sigma1);
					y33[tmpIndex] = y33Tilde * lambda / (lambda + 2 * sigma1);



					//_inputY3(tmpIndex) = std::max(-1.0f, std::min(1.0f, _inputY3(tmpIndex) + sigma2[tmpIndex] * (KX3(tmpIndex) + KX3(tmpIndex) - KX3old(tmpIndex) + ut[tmpIndex])));

					// where to get the u_t, 
					//y3[tmpIndex] = std::max(-1.0f, std::min(1.0f, y3[tmpIndex] + sigma2[tmpIndex] * (Kx3[tmpIndex] + Kx3[tmpIndex] - Kx3Old[tmpIndex] + ut[tmpIndex])));
					y4[tmpIndex] = std::max(-1.0f, std::min(1.0f, y4[tmpIndex] + sigma2[tmpIndex] * (Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old[tmpIndex] + ut[tmpIndex])));
					//y4[tmpIndex] = std::max(-1.0f, std::min(1.0f, y4[tmpIndex] + sigma2[tmpIndex] * (Kx4[tmpIndex] + Kx4[tmpIndex] - Kx4Old + ut[tmpIndex])));

					if (iterations % 50 == 0)
					{

						//d += std::abs((Y1old(tmpIndex) - _inputY1(tmpIndex)) / sigma1 - KX1old(tmpIndex) + KX1(tmpIndex)) +
						//	std::abs((Y1old(tmpIndex + nPx) - _inputY1(tmpIndex + nPx)) / sigma1 - KX1old(tmpIndex + nPx) + KX1(tmpIndex + nPx)) +
						//	std::abs((Y2old(tmpIndex) - _inputY2(tmpIndex)) / sigma1 - KX2old(tmpIndex) + KX2(tmpIndex)) +
						//	std::abs((Y2old(tmpIndex + nPx) - _inputY2(tmpIndex + nPx)) / sigma1 - KX2old(tmpIndex + nPx) + KX2(tmpIndex + nPx)) +
						//	std::abs((Y3old(tmpIndex) - _inputY3(tmpIndex)) / sigma2[tmpIndex] - KX3old(tmpIndex) + KX3(tmpIndex));

						d += std::abs((y11Old[tmpIndex] - y11[tmpIndex]) / sigma1 - Kx11Old[tmpIndex] + Kx11[tmpIndex]) +
							std::abs((y12Old[tmpIndex] - y12[tmpIndex]) / sigma1 - Kx12Old[tmpIndex] + Kx12[tmpIndex]) +
							std::abs((y13Old[tmpIndex] - y13[tmpIndex]) / sigma1 - Kx13Old[tmpIndex] + Kx13[tmpIndex]) +

							std::abs((y21Old[tmpIndex] - y21[tmpIndex]) / sigma1 - Kx21Old[tmpIndex] + Kx21[tmpIndex]) +
							std::abs((y22Old[tmpIndex] - y22[tmpIndex]) / sigma1 - Kx22Old[tmpIndex] + Kx22[tmpIndex]) +
							std::abs((y23Old[tmpIndex] - y23[tmpIndex]) / sigma1 - Kx23Old[tmpIndex] + Kx23[tmpIndex]) +

							std::abs((y31Old[tmpIndex] - y31[tmpIndex]) / sigma1 - Kx31Old[tmpIndex] + Kx31[tmpIndex]) +
							std::abs((y32Old[tmpIndex] - y32[tmpIndex]) / sigma1 - Kx32Old[tmpIndex] + Kx32[tmpIndex]) +
							std::abs((y33Old[tmpIndex] - y33[tmpIndex]) / sigma1 - Kx33Old[tmpIndex] + Kx33[tmpIndex]) +

							std::abs((y4Old[tmpIndex] - y4[tmpIndex]) / sigma2[tmpIndex] - Kx4Old[tmpIndex] + Kx4[tmpIndex]);


					}
				}

				if (iterations % 50 == 0)
				{
					err = (d*d + p * p) / (float)nPx;
				}

				if (iterations % 1000 == 0)
				{

					cout << "Iteration: " << iterations << " Residual " << err << endl;

					//mexPrintf("Iteration %d,Residual %e\n", iterations, err);
					//mexEvalString("drawnow;");
				}
			}
		}
		//cout << " in the L1TVOpticalFlowNonlinear 5 " << endl;

		//write output
		//#pragma omp parallel for
		//for (int tmpIndex = 0; tmpIndex < nPx; ++tmpIndex)
		//{

		// end of warps 

		//#pragma omp parallel for // collapse(3)
		//for (int k = 0; k < dim[2]; ++k)
		//{
		//	for (int j = 0; j < dim[1]; ++j)
		//	{
		//		for (int i = 0; i < dim[0]; ++i)
		//		{

#pragma omp parallel for
		for (int k = m_roi[4]; k < dim[2] - m_roi[5]; ++k)
		{
			for (int j = m_roi[2]; j < dim[1] - m_roi[3]; ++j)
			{
				for (int i = m_roi[0]; i < dim[0] - m_roi[1]; ++i)
				{
					int tmpIndex = index3DtoLinear(dim, i, j, k);
					

					_inputY1(i, j, k, 0) = (float)y11[tmpIndex];
					_inputY1(i, j, k, 1) = (float)y12[tmpIndex];
					_inputY1(i, j, k, 2) = (float)y13[tmpIndex];

					_inputY2(i, j, k, 0) = (float)y21[tmpIndex];
					_inputY2(i, j, k, 1) = (float)y22[tmpIndex];
					_inputY2(i, j, k, 2) = (float)y23[tmpIndex];

					_inputY3(i, j, k, 0) = (float)y31[tmpIndex];
					_inputY3(i, j, k, 1) = (float)y32[tmpIndex];
					_inputY3(i, j, k, 2) = (float)y33[tmpIndex];

					_inputY4(i, j, k, 0) = (float)y4[tmpIndex];

					//	_inputY(tmpIndex + 5 * nPx) = (double)y4[tmpIndex];
					//	_inputY(tmpIndex + 6 * nPx) = (double)y5[tmpIndex];

					//	_inputV[tmpIndex] += (double)v1[tmpIndex];
					//	_inputV[tmpIndex +  nPx] += (double)v2[tmpIndex];
					_inputV(i, j, k, 0) = (float)v1[tmpIndex];
					_inputV(i, j, k, 1) = (float)v2[tmpIndex];
					_inputV(i, j, k, 2) = (float)v3[tmpIndex];
					//	cout << _inputV(i, j, k, 0) << " " << _inputV(i, j, k, 1) << " " << _inputV(i, j, k, 2) << endl;

					/*		Outv1[tmpIndex] = _inputV(i, j, k, 0);
					Outv2[tmpIndex] = _inputV(i, j, k, 1);
					Outv3[tmpIndex] = _inputV(i, j, k, 2);*/
					//	_inputV(tmpIndex)
					//cout << Outv1[tmpIndex] << " ";
				}

			}
		}
		cout << " in the L1TVOpticalFlowNonlinear  5 " << endl;

		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cout << "=========Total time is :" << elapsed_secs << " s ======" << endl;
		int num1 = _inputV.width();
		int num2 = _inputV.height();
		int num3 = _inputV.depth();
		





		delete[] ux;
		delete[] uy;
		delete[] uz;
		delete[] ut;


		delete[] y11;
		delete[] y12;
		delete[] y13;
		delete[] y21;
		delete[] y22;
		delete[] y23;
		delete[] y31;
		delete[] y32;
		delete[] y33;
		delete[] y4;


		delete[] y11Old;
		delete[] y12Old;
		delete[] y13Old;
		delete[] y21Old;
		delete[] y22Old;
		delete[] y23Old;
		delete[] y31Old;
		delete[] y32Old;
		delete[] y33Old;


		delete[] tableI;
		delete[] tableJ;
		delete[] tableK;

		delete[] image1f;
		delete[] image2f;

		delete[] sigma2;
		//delete[] sigma3;
		//delete[] sigma4;
		delete[] tau1;
		delete[] tau2;
		delete[] tau3;

		delete[] v1;
		delete[] v2;
		delete[] v3;
		delete[] v1Old;
		delete[] v2Old;
		delete[] v3Old;



		delete[] y4Old;
		delete[] y5Old;

		delete[] Kty1;
		delete[] Kty2;
		delete[] Kty3;

		delete[] Kty1Old;
		delete[] Kty2Old;
		delete[] Kty3Old;

		delete[] Kx11;
		delete[] Kx12;
		delete[] Kx13;
		delete[] Kx21;
		delete[] Kx22;
		delete[] Kx23;
		delete[] Kx31;
		delete[] Kx32;
		delete[] Kx33;
		delete[] Kx4;
		/*delete[] Kx3;
		delete[] Kx4;
		delete[] Kx5;*/

		delete[] Kx11Old;
		delete[] Kx12Old;
		delete[] Kx13Old;
		delete[] Kx21Old;
		delete[] Kx22Old;
		delete[] Kx23Old;
		delete[] Kx31Old;
		delete[] Kx32Old;
		delete[] Kx33Old;
		//delete[] Kx4Old;
		delete[] Kx4Old;
		delete[] Kx5Old;

		//return ;
	}






	template< typename real, typename image >
	image L1TVOpticalFlowNonlinear_MultiScale3D(int *dim, const image &image1, const image &image2, image &_inputV,
		image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4, real _tol, real _lambda, int _maxIterations, int _norm, int _numberOfWarps,
		real eta, real sigma, int scales, real _huber, int *cut)
	{


	
		image I1_blur = image1.get_blur(sigma);
		image I2_blur = image2.get_blur(sigma);

	
		int size = image1.width()*image1.height()*image1.depth();
		
		imagelist Im1s;
		imagelist Im2s;
		imagelist inputVs;
		
		imagelist inputY1s;
		imagelist inputY2s;
		imagelist inputY3s;
		imagelist inputY4s;




		//imagelist dims;
		int twidth = I1_blur.width();
		int theight = I1_blur.height();
		int tdepth = I1_blur.depth();
		float t_eta = 1.0f;

		int tdims[] = { twidth ,theight, tdepth };
	
		for (int iscale = 0; iscale < scales; iscale++)
		{



			Im1s.insert(I1_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));
			Im2s.insert(I2_blur.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));
			inputVs.insert(_inputV.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY1s.insert(_inputY1.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY2s.insert(_inputY2.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY3s.insert(_inputY3.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 3, interpo));
			inputY4s.insert(_inputY4.get_resize(twidth*t_eta, theight*t_eta, tdepth*t_eta, 1, interpo));

			t_eta = t_eta * eta;


		}
		Im1s.reverse();
		Im2s.reverse();
		inputVs.reverse();
		//inputV2s.reverse();
		//inputV3s.reverse();
		//inputYs.reverse();
		inputY1s.reverse();
		inputY2s.reverse();
		inputY3s.reverse();
		inputY4s.reverse();
		//cout << "Test 1" << endl;
		////for (CImgList<>::iterator it = list.begin(); it<list.end(); ++it) (*it).mirror('x');
		// reture at last second scale;



		//	for (int iscale = 0; iscale < inputVs.size(); iscale++)
		for (int iscale = 0; iscale < inputVs.size(); iscale++)
		{
			int levels = inputVs.size() - iscale;
			cout << "At levels " << levels << endl;
			int _dim[] = { Im1s[iscale].width(),Im1s[iscale].height(),Im1s[iscale].depth() };
			//int t_dim[] = { inputY1s[iscale].width(),inputY1s[iscale].height(),inputY1s[iscale].depth() };
			cout << "_dim[]  :" << _dim[0] << " " << _dim[1] << " " << _dim[2] << endl;

			real _aabb[] = { 0,tdims[0], 0,tdims[1], 0,tdims[2] };

			//real _aabb[] = { -tdims[0]*0.5,tdims[0]*0.5, -tdims[1]*0.5,tdims[1]*0.5, -tdims[2]*0.5,tdims[2]*0.5 };
			cout << "tdims[]  :" << tdims[0] << " " << tdims[1] << " " << tdims[2] << endl;
			//if ((inputVs.size() - iscale) == 1)
			//{

			ROI[0] = (float)cut[0] * ((float)_dim[0] / (float)tdims[0]);
			ROI[1] = (float)cut[1] * ((float)_dim[0] / (float)tdims[0]);
			ROI[2] = (float)cut[2] * ((float)_dim[1] / (float)tdims[1]);
			ROI[3] = (float)cut[3] * ((float)_dim[1] / (float)tdims[1]);
			ROI[4] = (float)cut[4] * ((float)_dim[2] / (float)tdims[2]);
			ROI[5] = (float)cut[5] * ((float)_dim[2] / (float)tdims[2]);

			cout << "ROI: " << ROI[0] << " " << ROI[1] << " " << ROI[2] << " " << ROI[3] << " " << ROI[4] << " " << ROI[5] << " " << endl;
			//}
			//	L1TVOpticalFlowNonlinear3D(_dim, Im1s[iscale], Im2s[iscale], inputV1s[iscale], inputV2s[iscale], inputV3s[iscale],
			//	inputY1s[iscale], inputY2s[iscale],inputY3s[iscale], inputY4s[iscale], _tol, _lambda, _maxIterations, _norm,_numberOfWarps,_aabb);
			//image test = inputVs[iscale];
			L1TVOpticalFlowNonlinear3D(_dim, Im1s[iscale], Im2s[iscale], inputVs[iscale],
				inputY1s[iscale], inputY2s[iscale], inputY3s[iscale], inputY4s[iscale], _tol, _lambda, _maxIterations, _norm, _numberOfWarps, _aabb, _huber, ROI);
		

			if (iscale < inputVs.size() - 1)
			{

				int _dims[] = { inputVs[iscale + 1].width(), inputVs[iscale + 1].height(), inputVs[iscale + 1].depth() };
				//			int _olddims[] = { inputVs[iscale].width(), inputVs[iscale].height(), inputVs[iscale].depth() };

				inputVs[iscale + 1] = inputVs[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);

				inputY1s[iscale + 1] = inputY1s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
				inputY2s[iscale + 1] = inputY2s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
				inputY3s[iscale + 1] = inputY3s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 3, interpo);
				inputY4s[iscale + 1] = inputY4s[iscale].get_resize(_dims[0], _dims[1], _dims[2], 1, interpo);


			}


			//cout << inputVs.size()<<" "<< inputVs(i).width() << " " << inputVs(i).height() << endl;

		}
		////cout << Im1s.size() << endl;
		////system("pause");


		return inputVs[inputVs.size() - 1];
	}







}  // flow namespace


#endif




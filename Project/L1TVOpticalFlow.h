//
//#pragma once
//#include<iostream>
//using namespace std;
//#include "CImg.h"
//#include <map>
////#include"vec3.h"
//
////#include "advection.h"
//using namespace cimg_library;
//
//namespace cl = cimg_library;
//
//#define cimg_use_openmp
//
//typedef float				real;
////typedef geom::vec3<real>	real3;
////typedef vec2<real>	real2;
//typedef cl::CImg<real>		image;
//typedef cl::CImg<int>		table;
//typedef CImgList<real>  imagelist;
////typedef struct {
////	map< string, image >	point_data;
////	map< string, image >  cell_data;
////} frame_data;
////namespace flow_field {
//
////template< typename real, typename image >
//	void L1TVOpticalFlowNonlinear3D(int *dim, const image & image1, const image & image2, image & _inputV, image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4,
//		float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarpsreal, real *_aabb, float huber);
//	//void L1TVOpticalFlowNonlinear_MultiScale(const int *dim, const image & image1, const image & image2,  image & _inputV,  image & _inputY,
//	//	float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps,float eta, float sigma, int levels);
//
//	//template< typename real, typename image >
//	void L1TVOpticalFlowNonlinear_MultiScale3D( int *dim, const image &image1, const image &image2, image &_inputV,
//		image &_inputY1, image &_inputY2, image &_inputY3, image &_inputY4, real _tol, real _lambda, int _maxIterations, int _norm, int _numberOfWarps,
//		real eta, real sigma, int scales, real huber); //CPL2TVOpticalFlow
//	//template< typename real, typename image >
//	//void L1TVOpticalFlowNonlinear_MultiScale(const int *dim, const image &image1, const image &image2, image &_inputV,
//	//	image &_inputY1, image &_inputY2, image &_inputY3, float _tol, float _lambda, int _maxIterations, int _norm, int _numberOfWarps,
//	//	float eta, float sigma, int scales);
//
////}
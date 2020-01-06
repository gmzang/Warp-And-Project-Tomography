//#pragma once
#ifndef _FRAME_H_
#define _FRAME_H_


#include"CImg.h"

namespace cl = cimg_library;
typedef float real;
typedef cl::CImg<real> image;




// key frame item structure
class Frame_item {
public:
	Frame_item()
	{

	};
	Frame_item(int time, image img, int sp, int ep, int flow)
	{
		fTime = time;
		fVolume = img;
		//fIsActivated = active;
		fStartproj = sp;
		fEndproj = ep;
		fFlow = flow;
	};
	Frame_item(int volsize[3])
	{
		fTime = 0;
		fVolume = image(volsize[0], volsize[1], volsize[2], 1, 0.0f);
		//fIsActivated = false;
		fStartproj = -1;
		fEndproj = -1;
	};
	void SetVolSize(int volsize[3])
	{
		fsize[0] = volsize[0];
		fsize[1] = volsize[1];
		fsize[2] = volsize[2];
	}

	void GetVolSize(int *volsize)
	{
		volsize[0] = fsize[0];
		volsize[1] = fsize[1];
		volsize[2] = fsize[2];
	}
	//int fTime;
	float fTime;
	image fVolume;
	image fFlow;
	int fprojRange;
	//	int index;
	//bool fIsActivated;
	int fStartproj; // start projection for frame k, if not k=-1
	int fEndproj;// end projection for frame k, if not k=-1;
	int fsize[3];
};





// for projection item
typedef struct {
	//float angle; // each projection, assign unique angle while scanning
	//float time; // index for the # of proj 
//	image projection;
	//int frameid;
	int leftframe;
	int rightframe;
	float rElapsed;
	//int inframe[20];
	//void Inframe(int perframe){ // which frame does the projection belongs to.
	//	frameid = (int)(perframe / time);
	//}
}Proj_frame;

//// for frame item
//typedef struct {
//	float time;
//	image volume;
////	int index;
//	bool isActivated;
//	int startproj; // start projection for frame k, if not k=-1
//	int endproj;// end projection for frame k, if not k=-1;
//}Frame_item;
typedef struct {
	int volsize[3];
}size_item;


// for flow item
typedef struct {
	int index;
	image flow;
	int units;
}Flow_item;



#endif
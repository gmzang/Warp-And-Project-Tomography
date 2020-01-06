// VolumeData.h

#ifndef _VolumeData_h
#define _VolumeData_h

#include <stdio.h>


#include "CImg.h"
#include"Frames.h"
#include<vector>
////#include <map>
////#include"tv.h"
////
using namespace cimg_library;
////


class VolumeData
{
public:
	
	VolumeData();
	VolumeData(int d1,int d2, int d3, float sp1,float sp2,float sp3,float *data);
	virtual ~VolumeData();

	void Clean();

	

	void GetDimensions(int *dims) const;
	void SetDimensions( int *dims);

	void GetSpacings(float *spacings) const;
	void SetSpacings(const float *spacings) ;
	
	void SetProjectionMode(bool mode)
	{
		m_IsBackProjectMode=mode;

	}
	bool GetProjectionMode(){return m_IsBackProjectMode;}
	//void* GetData() const;
	//float * GetData(int iframe) const;

	//float * GetData(int iframe) const;

	CImg<float>  GetData(int iframe) ;
	//CImg<float> GetData();
	//float * GetData()const;
	//void SetData(float *_data) ;
	//void SetData(CImgList<float> _vol_list);
	void SetData();
	void setParas(int *vols,double vs, int nframes, int totalproj);
	

	
	//void WriteStructuredFile(const char* filename) const;
	void ReadStructuredFile();
	//void ReadStructuredFile(const char* filename);
	void _allocate();
	bool m_IsBackProjectMode;
	//float *m_EachVoxel_TotalWeight;
	//float *m_Voxel_Contributions;
	//float *m_correction;

	int m_Dims[3];
	float m_Spacings[3];

	//void *m_Data;
	//float *m_Data;
	int m_volX;
	int m_volY;
	int m_volZ;
	double Vol_Spacing;
	
	std::vector<Frame_item> m_Frameslist[5];
	int m_fInLevel[5];
	float m_tUnitInLevel[5];

	std::vector<size_item> m_Volumesize;
	std::vector<Proj_frame> m_ProjIndex;
	int m_nframes;
	int m_totalProjs;
	float m_dowmsampleFactor;
	int m_roi[8];
	int m_levels;
	int m_base;
};

#endif

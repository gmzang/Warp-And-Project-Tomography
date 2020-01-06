//PSART.h

#ifndef __PSART_h
#define __PSART_h

//#define cimg_use_openmp
#define cimg_use_tiff
#define cimg_use_tif
#include "CImg.h"
#include "Scene.h"
#include <string>
#include"Frames.h"
class VolumeData;
class RayCaster;
class TypicalRayProducer;
class CoherentVolume;
//class GridPartition;



class PSART : public Scene
{
public:

	

	PSART(double d1, double d2, double d3, const char *outputfile, int *imgwh, int projNo, int iters, int sartiters,
		double alp, double sid, const char **projsfile, double *oxyz,  double sdd, double ds, const char * _prior, 
		const char * _bp, int *volxyz,int nframes,double *startdg, int nrounds);
	//PSART(double d1,double d2,double d3, const char *outputfile,int imgw,int imgh,int projNo,int iters,int sartiters, 
	//	double alp,double sid,const char *projsfile,double ox, double oy, double oz,double sdd,double offx,double offy,double ds,const char * _prior,const char * _bp, int *test, int *test2);
	virtual ~PSART();
		void InitViews();
	bool RotateVolume(int i );
	// data access 
	void SetData(VolumeData *data);
	void SART(int iframes, int level, int *dimension, int *projsize, int nb_sart);//AdaptiveWarp
	void SART_Proximal_Operator(int iframes,  int nb_sart, int outloop,int stIter, int flevel);
	void AdaptiveWarp(int iframes, int level, int *dimension, int *projsize, int nb_sart, int outloop, float *aabb);

	//void ShrinkProjs(int out_loop);
	void WarpAndProjct(int iproj, int iframes, int *dimension, float rtime,  float *aabb, int outloop, int warpiter,  char * _file_header, int flevel);
	void SetBackGroundColor(float r,float g, float b);
	void GetBackGroundColor(float& r,float& g,float& b);

	void SetIsovalue(float isovalue);
	float GetIsovalue();

	void SetSurfaceColor(float r,float g,float b);
	void GetSurfaceColor(float &r,float &g,float &b);

	void SetLightDirection(float x,float y,float z);
	void SetLightIntensity(float intensity);
	void SetLightColor(float r,float g,float b);

	void SetAmbient(float value);
	float GetAmbient();

	void SetDiffuse(float value);
	float GetDiffuse();

	void SetSpecular(float value);
	float GetSpecular();

	void SetSpecularPower(float value);
	float GetSpecularPower();
	
	void  SetSampleDistance(float sd);
	float GetSampleDistance();

	
	virtual void Init();
	virtual void Resize(int width, int height);
	//void addMhaHeader(const int volDims[3],const float volSpacings[3]);
protected:
	
	//virtual void RenderContent(int iframe,int level);
	virtual void RenderContent(int  iframe, bool isWarped, int ProIdx, float rtime, int flevel);
	virtual void GetSize3D(float size[3]);
	virtual void Render();
	
private:
	// main loop
	//bool _renderCoreFuc(int iframe, int level);
	bool _renderCoreFuc(int iframe, bool isWarpped, int ProIdx, float rtime, int flevel);
	//bool _renderCoreFuc();

	// data
	VolumeData* m_Data;
	

	// ray caster
	RayCaster* m_RayCaster;
	// ray generator
	TypicalRayProducer* m_RayProducer;
	// main direction
	int m_mainDir;
	int m_projNo;
	float m_BackgroundColor[3];
	//float *VoxelValue;
	float **pImg;

	float *pProImg;

	Quaternion *pQuatProjs;
	int m_nframes;
	int m_nrounds;

	//float ** pY;
	//float *pXbar;
	//float *pXPrev;
	double m_s;
	double m_t;
	double m_l;
	const char *pOutput;
	//const char *pInput;
	 const char **pInput;
	//float *pD;
	//float *pZ_SAD;
	double m_mu;
	float **m_startDegree;
	float m_huberFactor;
	double m_volspace;
	int m_voldims[3];


	
//		GridPartition *m_GridPartition;
	//CoherentVolume* m_Volume;
};

#endif
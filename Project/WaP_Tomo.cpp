//PSART.cpp

#include "WaP_Tomo.h"
#include "RayCaster.h"
#include "TypicalRayProducer.h"
#include "VolumeData.h"
#include <iostream>
#include <malloc.h>
#include <stdlib.h>
#include <string>
#include "GL/glew.h"
#include "ImageIO.h"
#include <malloc.h>
#include <fstream>
#include <time.h>
#include <math.h>
#include <sstream>
#include "tiffio.h"
#include "tiff.h"

#include"optical_flow.h"
#include"advection.h"



#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <string>
#include<Eigen/Core>
#include <Eigen/Dense>
#include<Eigen/SVD>
#include <Eigen/Householder>
#include <time.h> 

#include "CoherentVolume.h"
//#include "GridPartition.h"
#include<fstream>
//#include "CImg.h"
#include <map>
#include"advection.h"
using namespace cimg_library;

namespace cl = cimg_library;

typedef cl::CImg<float>		image;
typedef cl::CImg<int>		table;
typedef CImgList<float>  imagelist;

image tim;
image imROI;



using namespace std;
int myrandom(int i) { return std::rand() % i; }
#if !defined(Circular_Angle)
# define Circular_Angle 360//or150 198.787878f 可能角度需要改变
//# define Angle_Step Circular_Angle/(PRO_NO-1)
#endif

#if !defined(IS_COUNTERCLOCK)
#  define IS_COUNTERCLOCK 1 //counterclockwise 。 1 for tif  1 change to 0 to see the rabbit  
// counterclockwise means from x-ray source, the object is scanned clockwise from top view 
#define STV   1
#define SAD   2
#define ATV   3
#define Voxelbased 1
#define Raybased  2
#endif

int PRO_WIDTH;
int PRO_HEIGHT;
int iminY;
static int PRO_NO = 6; //
int scannerRound = 32;  // change place 2
//static int PRO_NO = 360; //
static int ITERATIONS = 1;//
static int CP_ITERS = 1;//
Eigen::MatrixXf pY_warp;
Eigen::MatrixXf pY;
Eigen::MatrixXf pY_Adaptive;

double m_theta = 1.0;
static const float pi = 3.1415926536f;

float m_para;
float G_Max = 10.0f;
double Alpha = 0.1;
double m_alpha = Alpha;
#define USING_PRIOR ATV;
int  Prior = ATV;
int  Backprojector = Raybased;
int  Joint_iter = 20;
//#  define Factor sqrtf(2*m_para) //change  //finish chang 1


float  L = 8.0;

float CP_tau = 1; //0.3
////
//////float ATV_sigma = .9 / (L*ATV_tau);
//float ATV_sigma = 0.3;
float CP_sigma = .9 / (L*CP_tau);
float	m_aabb[6];

float TV_w = 0.01; // 
float TP_w = 0.05;  // 
float BC_w = 0.1;  //
float SART_lambda = CP_tau / 2;
#  define Factor sqrtf(2*SART_lambda) 

int cut[] = { 8,8,8,8,8,8,8,8 };
int boundx1 = 5;
int boundx2 = 5;
int boundy1 = 5;
int boundy2 = 5;
static inline void _Cal_Norm(float &u, float sigma)
{
	/*if(u>m_sigma) u=m_sigma;
	else if(u<-m_sigma) u=-m_sigma;*/
	if (u > sigma) u = sigma;
	else if (u < -sigma) u = -sigma;

}
//int dims[3];
int timestorender = 0;
int rot = 1;
//static int _index_= 0;
static int _ReadyToRotate_ = 1;

static double Angle_Step = (double)Circular_Angle / (double)(PRO_NO); //test if identical to RTK results


static int projsEachRound = 6;
//static int projsEachRound = 1;
static double degreeEachProj = 360.0 / projsEachRound;
double inputdegree[3600];// number of projection


float bak_densityaverage = 0.0f;
bool SetOriginalPos = false;
int index_of_project = 0;
char* pOutputFile;
const char* outputfilename;
int len;
string output = "t";  //RTK_Origin_0  SL128_RTK_KB
string outputformat = ".mha";
int numofSART = 0;
//int NO_SART=1;
int *idxSART;
int *idxSART_Proximal;
int *idxAdaptive;
int noindex = 0;
float Views = 0;
double sampledistance = 0.0f;
int vh_x1 = 0;
int vh_x2 = 0;
int vh_y1 = 0;  //120   20
int vh_y2 = 0;
int vh_z1 = 0;
int vh_z2 = 0;
int CP_iter = 0;
image warpI2;
int idx_proj = 0;



#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#ifndef idx
#define idx(i,j,twidth)            (i + j*twidth)
#endif

#ifndef idx3D
#define idx3D(i,j,k,twidth,theight)            (i + j*twidth+ k*twidth*theight)
#endif


#if !defined(TestIter)
#  define TestIter  63 //45
#endif


int levels = 1;
imagelist *projslevel = new imagelist[levels];





template<typename _Matrix_Type_>
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon())
{
	Eigen::JacobiSVD< _Matrix_Type_ > svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}


float inline cubicInterpolate(float p[4], float x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x * (2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x * (3.0*(p[1] - p[2]) + p[3] - p[0])));

}


float inline bicubicInterpolate(float p[4][4], float x, float y) {
	float arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);

}





PSART::PSART(double d1, double d2, double d3, const char *outputfile, int *imgwh, int projno, int iters, int sartiters,
	double alp, double sid, const char **projsfile, double *oxyz, double sdd, double ds, const char * _prior, const char * _bp, int *volxyz,
	int nframes, double *startdg, int nrounds)
{

	m_nframes = nframes;
	m_nrounds = nrounds;

	m_RayProducer = new TypicalRayProducer(volxyz[0], volxyz[1], volxyz[2], imgwh[0], imgwh[1], ds, _bp);
	Alpha = alp;
	m_alpha = Alpha;

	PRO_WIDTH = imgwh[0];
	PRO_HEIGHT = imgwh[1];
	iminY = imgwh[0];
	
	TV_w = d1;
	TP_w = d2;
	BC_w = 0.0f;
	if (strcmp(_bp, "Raybased") == 0)
	{
		Backprojector = Raybased;
	}
	else
		Backprojector = Voxelbased;


	m_halfFrustumWidth = (imgwh[0] - 1)*ds*0.5f;
	m_halfFrustumHeight = (imgwh[1] - 1)*ds*0.5f;

	m_projNo = projno;
	Angle_Step = (double)Circular_Angle / (double)(PRO_NO);
	Views = Angle_Step;




	ITERATIONS = sartiters;
	CP_ITERS = iters;

	cout << "angles:  " << degreeEachProj << " " << Angle_Step << endl;

	m_sid = sid;
	m_sdd = sdd;
	m_ds = ds;
	
	m_ViewPort[0] = imgwh[0];
	m_ViewPort[1] = imgwh[1];



	m_Translate[0] = -oxyz[0];
	m_Translate[2] = +oxyz[1];
	//m_Translate[2] = -oxyz[1];
	m_Translate[1] = -oxyz[2];
	pOutput = outputfile;

	//pOutput=projs;
	//	pInput=projsfile;


	//pQuatProjs = new Quaternion[projno];

	m_Data = NULL;
	m_BackgroundColor[0] = m_BackgroundColor[1] = m_BackgroundColor[2] = 0.0f;

	m_RayCaster = new RayCaster;

	pInput = new const char*[5000];



	idxSART = new int[5000];
	idxSART_Proximal = new int[5000];
	idxAdaptive = new int[5000];

	cout << "Test ---2" << endl;
	//ifstream myfile("lds.txt");

	int idxrounds = 0;
	//for (int i = 0;i < nframes;i++)
	for (int i = 0;i < 1000;i++)
	{

	
		pInput[i] = new char[PRO_NO];
		// need check this initilization
		pInput[i] = projsfile[i];
		//pQuaternion[i] = new Quaternion[PRO_NO];
		idxSART[i] = 1;
		idxSART_Proximal[i] = 1;
		idxAdaptive[i] = 1;
	}
	cout << "Test ---3" << endl;


	ifstream myinput("ours_foam_compress.txt");
	//
	for (int i = 0;i < m_projNo;i++)
	{
		myinput >>  inputdegree[i];
		//cout << i<<" "<< inputdegree[i] << " ,";
	}
	myinput.close();


	// test all parameters here:
	cout << PRO_WIDTH << " " << PRO_HEIGHT << " " << volxyz[0] << " " << volxyz[1] << " " << volxyz[2] << endl;

	

	projslevel[0].insert(image(pInput[0]).get_split('z', m_projNo));


	float t_min = 10000.0f;
	float t_max = 0.0f;
	string t_str = pInput[0];
	cout << t_str << endl;
	//cout << "input volume max:	" << _max << "	_min:" << _min << endl;
	t_max = projslevel[0][0].max_min(t_min);
	G_Max = t_max/1.0f;
	cout << "G_Max:	" << G_Max << endl;



}

PSART::~PSART()
{
	delete m_RayProducer;
	delete m_RayCaster;
}

void PSART::SetData(VolumeData *data)
{
	//	cout << "testing 2" << endl;
	m_Data = data;
	SetSize3DModified();

	if (!data) return;

	float spacings[3];
	//m_Data->GetDimensions(dims);
	m_Data->GetDimensions(m_voldims);
	m_Data->GetSpacings(spacings);
	//	new added  GUANGMING 10-24 Delete it.
	//m_Volume->SetData(m_Data);
//	m_GridPartition->Init(m_Volume,0.4f);

	m_volspace = spacings[1];
	m_RayProducer->SetVolSize(m_voldims, spacings[1]);

	m_RayCaster->SetSampleDistance(sqrt((spacings[0]*spacings[0]+spacings[1]*spacings[1]+spacings[2]*spacings[2]))/2);
	// change the sample rate
	//m_RayCaster->SetSampleDistance(sqrt(spacings[0] * spacings[0] + spacings[1] * spacings[1] + spacings[2] * spacings[2]));
	sampledistance = m_RayCaster->GetSampleDistance();


	cout << "Writing output file dims: " << m_voldims[0] << " " << m_voldims[1] << " " << m_voldims[2] << " " << spacings[0] << " " << spacings[1] << " " << spacings[2] << endl;




	imROI = image(m_voldims[0], m_voldims[1], m_voldims[2], 1.0);
#pragma omp parallel for
	for (int k = 0; k < m_voldims[2]; k++)
		for (int j = 0; j < m_voldims[1]; j++)
			for (int i = 0; i < m_voldims[0]; i++)
			{
				{
					{

						imROI(i, j, k) = 1.0f;
						if (((i < m_Data->m_roi[6]) && (k < m_Data->m_roi[7])) || ((i < m_Data->m_roi[6]) && (k > (m_voldims[2] - m_Data->m_roi[7]))))
							imROI(i, j, k) = 0.0f;
						if (((i > (m_voldims[0] - m_Data->m_roi[6])) && (k < m_Data->m_roi[7])) || ((i > (m_voldims[0] - m_Data->m_roi[6])) && (k > (m_voldims[2] - m_Data->m_roi[7]))))
							imROI(i, j, k) = 0.0f;

						if (((k < m_Data->m_roi[6]) && (i < m_Data->m_roi[7])) || ((k < m_Data->m_roi[6]) && (i > (m_voldims[0] - m_Data->m_roi[7]))))
							imROI(i, j, k) = 0.0f;
						if (((k > (m_voldims[2] - m_Data->m_roi[6])) && (i < m_Data->m_roi[7])) || ((k > (m_voldims[2] - m_Data->m_roi[6])) && (i > (m_voldims[0] - m_Data->m_roi[7]))))
							imROI(i, j, k) = 0.0f;

						if (i<m_Data->m_roi[0] || i>m_voldims[0] - m_Data->m_roi[1] || j<m_Data->m_roi[2] || j>m_voldims[1] - m_Data->m_roi[3] || k<m_Data->m_roi[4] || k>m_voldims[2] - m_Data->m_roi[5])
							imROI(i, j, k) = 0.0f;

						//else 
					}
				}
			}

	m_aabb[0] = 0;
	m_aabb[1] = m_voldims[0];
	m_aabb[2] = 0;
	m_aabb[3] = m_voldims[1];
	m_aabb[4] = 0;
	m_aabb[5] = m_voldims[2];


	int psize = m_voldims[0] * m_voldims[1] * m_voldims[2];




	cout << "Finish InitViews()" << endl;

}

void PSART::SetBackGroundColor(float r, float g, float b)
{
	m_BackgroundColor[0] = r; m_BackgroundColor[1] = g; m_BackgroundColor[2] = b;
}

void PSART::GetBackGroundColor(float& r, float& g, float& b)
{
	r = m_BackgroundColor[0]; g = m_BackgroundColor[1]; b = m_BackgroundColor[2];
}



void PSART::Init()
{
	glewInit(); // initial glew
}

void PSART::Resize(int width, int height)
{
	Scene::Resize(width, height);
	cout << "width: " << width << " height: " << height << endl;
	m_RayProducer->SetViewPort(width, height);
}

void  PSART::SetSampleDistance(float sd)
{
	m_RayCaster->SetSampleDistance(sd);
}
float PSART::GetSampleDistance()
{
	return m_RayCaster->GetSampleDistance();
}





void PSART::Render()
{
	

	cout << "In Render()" << endl;
	int m_outloop[] = { 1,1,1,1,1 };// { 2, 5, 7, 9 };//5
	//int m_outloop[] = { 2,3,3,4,4 };// 
	int nb_init_sart = 15;
	int timeofPO = 1;
	int m_each_level_iter = 3; //3
	image _warpI1;
	image _warp_iframe_plus;
	int nb_proxiSART = 2 ; // 6
	int nb_adaptiveWarp = 1;
	int nb_projwarp = 1; //4
	int test_threshold = 5;

	int nb_volOpt[] = { 2,10,21,6,21,15 }; //bubble

	int nlevels = 2;

	int dim[] = { m_voldims[0], m_voldims[1], m_voldims[2] };
	float aabb[] = { 0,dim[0], 0,dim[1], 0,dim[2] };

	for (int ilistLevel = 0; ilistLevel < nlevels; ilistLevel++)
	{

	m_nframes = m_Data->m_fInLevel[ilistLevel];
	cout << "Current m_nframes: " << m_nframes << endl;

	
	imagelist Y_tv(m_nframes + 3, m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);//(dims[0], dims[1], dims[2], 3, 0.0);
				   // remove 6
	imagelist Y_tp(m_nframes + 3, m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);// temporal prior//(dims[0], dims[1], dims[2], 3, 0.0);
	imagelist Y_bc(m_nframes + 3, m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);// temporal prior//(dims[0], dims[1], dims[2], 3, 0.0);

	imagelist X_bar(m_nframes + 3, m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);//(dims[0], dims[1], dims[2], 1, 0);
	imagelist X_prev(m_nframes + 3, m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);//(dims[0], dims[1], dims[2], 1, 0);
	// can change to curcurrent function
	//cout << " size: " << X_bar.size() << " " << X_bar.size() << " " << Y_bc.size() << endl;
	//system("pause");
	for (int outloop = 0; outloop < m_outloop[ilistLevel]; outloop++)
	{
		cout << "outloop:	" << outloop<< endl;

	
		// Initialization
		image tim;
		image imROI;
		image tmpimage;
		m_huberFactor = 1.0f;
		float dt = 1.0;


		idx_proj = 0;

		m_Data->m_roi[0] = (float)cut[0] * ((float)m_voldims[0] / (float)m_Data->m_Volumesize[0].volsize[0]);
		m_Data->m_roi[1] = (float)cut[1] * ((float)m_voldims[0] / (float)m_Data->m_Volumesize[0].volsize[0]);
		m_Data->m_roi[2] = (float)cut[2] * ((float)m_voldims[1] / (float)m_Data->m_Volumesize[0].volsize[1]);
		m_Data->m_roi[3] = (float)cut[3] * ((float)m_voldims[1] / (float)m_Data->m_Volumesize[0].volsize[1]);
		m_Data->m_roi[4] = (float)cut[4] * ((float)m_voldims[2] / (float)m_Data->m_Volumesize[0].volsize[2]);
		m_Data->m_roi[5] = (float)cut[5] * ((float)m_voldims[2] / (float)m_Data->m_Volumesize[0].volsize[2]);
		m_Data->m_roi[6] = (float)cut[6] * ((float)m_voldims[2] / (float)m_Data->m_Volumesize[0].volsize[2]);
		m_Data->m_roi[7] = (float)cut[7] * ((float)m_voldims[2] / (float)m_Data->m_Volumesize[0].volsize[2]);

		





		// if the initialization is correct.
		pY_warp = Eigen::MatrixXf::Zero(1, PRO_WIDTH*PRO_HEIGHT);
		pY_Adaptive = Eigen::MatrixXf::Zero(1, PRO_WIDTH*PRO_HEIGHT);
		pY = Eigen::MatrixXf::Zero(PRO_NO, PRO_WIDTH*PRO_HEIGHT);


		cout << "==============================================" << endl;
		cout << "S1: Frames-Reconstruction-subOptimization starts " << endl;
		cout << "==============================================" << endl;
		// step 1: Frames-Reconstruction 
		cout << "======= weight for TV: " << TV_w << ", and TP: " << TP_w << endl;


		char basename[1024];
		//char flow_field[1024];
		//sprintf(out_dir, ".");
		sprintf(basename, "FirstStep");


		for (int i_cp_vol = 1;i_cp_vol <= nb_volOpt[ilistLevel];i_cp_vol++)
		{
			cout << "=============" << timeofPO++ << "  times in Proximal Operation iteration ==============" << endl;
			int kkk = 0;
			int inY = m_voldims[0];
			int inZ = m_voldims[0] * m_voldims[1];
			int inX = 1;

			for (int iframe = 0; iframe < m_nframes; iframe++)
			{
				cout << "iframe: " << iframe << endl;
		
				// Dual part 
				if ( (i_cp_vol == 0) && outloop==0)
				{
					X_bar[iframe] = m_Data->m_Frameslist[ilistLevel][iframe].fVolume;
					X_bar[iframe + 1] = m_Data->m_Frameslist[ilistLevel][iframe + 1].fVolume;
					//X_bar[iframe] = m_Data->m_Volumelist[iframe];
					//X_bar[iframe + 1] = m_Data->m_Volumelist[iframe + 1];
				}
				cout << m_voldims[0] << " " << m_voldims[1] << " " << m_voldims[2] <<" "<< m_Data->m_Frameslist[ilistLevel].size()<< endl;
				//_warp_iframe_plus = advection::advect(m_voldims, aabb, m_Data->m_Frameslist[ilistLevel][iframe].fFlow, X_bar[iframe + 1], -dt);

				Y_tp[iframe] += CP_sigma * (X_bar[iframe + 1] - X_bar[iframe]);

				Y_bc[iframe] += CP_sigma * (_warp_iframe_plus - X_bar[iframe]);

				#pragma omp parallel for
				for (int k = m_Data->m_roi[4]; k < m_voldims[2] - m_Data->m_roi[5]; k++) {
					for (int j = m_Data->m_roi[2]; j < m_voldims[1] - m_Data->m_roi[3]; j++) {
						for (int i = m_Data->m_roi[0]; i < m_voldims[0] - m_Data->m_roi[1]; i++) {

							/*	if (i<cut[0] || i>m_voldims[0] - cut[1] || j<cut[2] || j>m_voldims[1] - cut[3] || k<cut[4] || k>m_voldims[2] - cut[5])
								{
									continue;
								}*/
							Y_tv[iframe](i, j, k, 0) += CP_sigma * (X_bar[iframe](i + 1, j, k) - X_bar[iframe](i, j, k))*m_huberFactor;
							Y_tv[iframe](i, j, k, 1) += CP_sigma * (X_bar[iframe](i, j + 1, k) - X_bar[iframe](i, j, k))*m_huberFactor;
							Y_tv[iframe](i, j, k, 2) += CP_sigma * (X_bar[iframe](i, j, k + 1) - X_bar[iframe](i, j, k))*m_huberFactor;

							Y_tv[iframe](i, j, k, 0) = max(-TV_w, min(TV_w, Y_tv[iframe](i, j, k, 0)));
							Y_tv[iframe](i, j, k, 1) = max(-TV_w, min(TV_w, Y_tv[iframe](i, j, k, 1)));
							Y_tv[iframe](i, j, k, 2) = max(-TV_w, min(TV_w, Y_tv[iframe](i, j, k, 2)));

							// remove 2
							//P164 , for L2 norm in dual part,equal to y=y*k3/(k3+2*ATV_sigma)
							Y_tp[iframe](i, j, k) = Y_tp[iframe](i, j, k)*TP_w / (TP_w + 2 * CP_sigma);

							//Y_tp[iframe](i, j, k)= max(-TP_w, min(TP_w, Y_tp[iframe](i, j, k)));
							Y_bc[iframe](i, j, k) = max(-BC_w, min(BC_w, Y_bc[iframe](i, j, k)));

							//Y_tp[iframe](i, j, k, 1) = max(-TP_w, min(TP_w, Y_tp[iframe](i, j, k, 1)));
							//Y_tp[iframe](i, j, k, 2) = max(-TP_w, min(TP_w, Y_tp[iframe](i, j, k,2)));

						}
					}
				}

				//X_prev[iframe] = m_Data->m_Volumelist[iframe];
				X_prev[iframe] = m_Data->m_Frameslist[ilistLevel][iframe].fVolume;


				// Primal part
				if (iframe == 0)
				{

					m_Data->m_Frameslist[ilistLevel][iframe].fVolume -= CP_tau * (-Y_tp[iframe]);
					m_Data->m_Frameslist[ilistLevel][iframe].fVolume -= CP_tau * (-Y_bc[iframe]);
				
				}
				else
				{
				
					m_Data->m_Frameslist[ilistLevel][iframe].fVolume -= CP_tau * (Y_tp[iframe - 1] - Y_tp[iframe]);
					m_Data->m_Frameslist[ilistLevel][iframe].fVolume -= CP_tau * (_warpI1 - Y_bc[iframe]);

				}



				#pragma omp parallel for
				for (int k = m_Data->m_roi[4]; k < m_voldims[2] - m_Data->m_roi[5]; k++) {
					for (int j = m_Data->m_roi[2]; j < m_voldims[1] - m_Data->m_roi[3]; j++) {
						for (int i = m_Data->m_roi[0]; i < m_voldims[0] - m_Data->m_roi[1]; i++) {

						

							float y1 = -(Y_tv[iframe](i, j, k, 0) - Y_tv[iframe](i - 1, j, k, 0) + Y_tv[iframe](i, j, k, 1) -
								Y_tv[iframe](i, j - 1, k, 1) + Y_tv[iframe](i, j, k, 2) - Y_tv[iframe](i, j, k - 1, 2));
							
							m_Data->m_Frameslist[ilistLevel][iframe].fVolume(i, j, k) -= (CP_tau*y1);;
							//m_Data->m_Data[idx(i, j, k)] -= (ATV_lambda*(Y_tv(i, j, k, 0) - Y_tv(i + 1, j, k, 0) + Y_tv(i, j, k, 1) - Y_tv(i, j + 1, k, 1) + Y_tv(i, j, k, 2) - Y_tv(i, j, k + 1, 2)));

							//m_Data->m_Data[p] = max(min(G_Max, m_Data->m_Data[p]), 0);
							m_Data->m_Frameslist[ilistLevel][iframe].fVolume(i, j, k) = max(min(G_Max, m_Data->m_Frameslist[ilistLevel][iframe].fVolume(i, j, k)), 0);
						}
					}
				}

				_warp_iframe_plus = 0.0f;
				//_warp_iframe_minus = 0.0f;
				_warpI1 = 0.0f;

				// change _5_
				SART_Proximal_Operator(iframe, nb_proxiSART, outloop, i_cp_vol, ilistLevel);

				
	
				idxSART_Proximal[iframe] = 1;
				idxAdaptive[iframe] = 1;


				#pragma omp parallel for
				for (int z = m_Data->m_roi[4];z < m_voldims[2] - m_Data->m_roi[5];z++)
				{
					for (int y = m_Data->m_roi[2];y < m_voldims[1] - m_Data->m_roi[3];y++)
					{
						for (int x = m_Data->m_roi[0];x < m_voldims[0] - m_Data->m_roi[1];x++)
						{
							
							X_bar[iframe](x, y, z) = 2 * m_Data->m_Frameslist[ilistLevel][iframe].fVolume(x, y, z) - X_prev[iframe](x, y, z);

						}

					}

				}
				

				cout << "=========Image reconstruction=====-" << CP_iter << "-CP iter " << iframe << "frame completed" << endl;
				//cout << endl;
			}//	for (int iframe = 0; iframe < m_nframes; iframe++)

			//if (outloop == 0)
			//	ShrinkProjs(outloop);
		}  //for (CP_iter = 0;CP_iter < CP_ITERS;CP_iter++)
	
		   
		   
		   




//cout << "==============================================" << endl;
//cout << "S2: FlowField-Estimation-subOptimization starts " << endl;
//cout << "==============================================" << endl;

float eta = 0.5;  //downsample eta 0.65 and sigma 0.5
float sigma = 0.8;  //blur 0.5
float precision = 0.008; //0.001 IS GOOD
image Y1 = image(m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);
image Y2 = image(m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);
image Y3 = image(m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);
image Y4 = image(m_voldims[0], m_voldims[1], m_voldims[2], 1, 0.0);
image uvw = image(m_voldims[0], m_voldims[1], m_voldims[2], 3, 0.0);
image I1;
image I2;
//float  smoothness = 0.5;
float  smoothness = 0.05; //0.03
int _norm = 1; //  huber==4
float huber = 1.0f;
int iters = 3000;
int scales = 5;
int nwarps = 3;


	//#pragma omp parallel for
	//for (int flowindx = 0; flowindx < m_nframes - 1; flowindx++)
	for (int flowindx = 0; flowindx < m_nframes - 1; flowindx++)
	{
		
		// reset values to 0;
		I1.fill(0.0f);
		I2.fill(0.0f);
		uvw.fill(0.0f);
		Y1.fill(0.0f);
		Y2.fill(0.0f);
		Y3.fill(0.0f);
		Y4.fill(0.0f);

	
		//I1 = m_Data->m_Frameslist[flowindx].fVolume;
		//I1 = advection::advect(dim, aabb, m_Data->m_Frameslist[flowindx].fFlow, m_Data->m_Frameslist[flowindx].fVolume, +dt);
		I1 = m_Data->m_Frameslist[ilistLevel][flowindx].fVolume;
		I2 = m_Data->m_Frameslist[ilistLevel][flowindx + 1].fVolume;
		uvw = m_Data->m_Frameslist[ilistLevel][flowindx].fFlow;
		

		m_Data->m_Frameslist[ilistLevel][flowindx].fFlow =optical_flow::L1TVOpticalFlowNonlinear_MultiScale3D(m_voldims, I1, I2, uvw,
			Y1, Y2, Y3, Y4, precision, smoothness, iters, _norm, nwarps, eta, sigma, scales, huber, m_Data->m_roi);
		
		char out_dir[1024];
		char basename[1024];
		char flow_field[1024];
		//sprintf(out_dir, ".");
		sprintf(basename, "VelocityField");
		//sprintf(flow_field, "flow_field");


		image _warpedI1 = advection::advect(m_voldims, aabb, m_Data->m_Frameslist[ilistLevel][flowindx].fFlow, I2, -dt);
		image _warpedI2 = advection::advect(m_voldims, aabb, m_Data->m_Frameslist[ilistLevel][flowindx].fFlow, I1, +dt);
		//image _warpedI2 = advection::advect(m_voldims, aabb, m_Data->m_Frameslist[flowindx].fFlow, m_Data->m_Frameslist[flowindx].fVolume, +dt);
		_warpedI1.save_analyze(str_format("warpedVolume/warpedI%03d.s%1.3f.h%3.3f.f%03d.outloop%03d.f%03d_back.hdr" , flowindx, smoothness, huber, flowindx,outloop, flowindx).c_str());
		_warpedI2.save_analyze(str_format("warpedVolume/warpedI%03d.s%1.3f.h%3.3f.f%03d.outloop%03d.f%03d_forward.hdr", flowindx + 1, smoothness, huber, flowindx, outloop, flowindx).c_str());

		m_Data->m_Frameslist[ilistLevel][flowindx].fFlow.get_channel(0).save_analyze(str_format("%s/Flow%03d.s%2.3f.h%1.3f.f%03d.outloop%03d.u.hdr", basename, dim[0], smoothness, huber, flowindx, outloop).c_str());
		m_Data->m_Frameslist[ilistLevel][flowindx].fFlow.get_channel(1).save_analyze(str_format("%s/Flow%03d.s%2.3f.h%1.3f.f%03d.outloop%03d.v.hdr", basename, dim[0], smoothness, huber, flowindx, outloop).c_str());
		m_Data->m_Frameslist[ilistLevel][flowindx].fFlow.get_channel(2).save_analyze(str_format("%s/Flow%03d.s%2.3f.h%1.3f.f%03d.outloop%03d.w.hdr", basename, dim[0], smoothness, huber, flowindx, outloop).c_str());
		//}
		//m_Data->m_Frameslist[ilistLevel][2].fFlow.get_acos()



	}



			/*
			Step 3: warp and project
			*/
		//cout << "==============================================" << endl;
		//cout << "Warp-And-Project-subOptimization starts " << endl;
		//cout << "==============================================" << endl;
		//// initilization


		char _Base[1024];
		//char flow_field[1024];
		//sprintf(out_dir, ".");
		sprintf(_Base, "WarpStep");

		
		for (int wpIter = 0; wpIter < nb_projwarp; wpIter++)
		{
			// We can use openmp
			//#pragma omp parallel for
			float trange = m_Data->m_tUnitInLevel[ilistLevel];
			for (int fidx = 0; fidx < m_Data->m_nframes; fidx++)
			{
				int tleft = m_Data->m_Frameslist[ilistLevel][fidx].fStartproj;
				int tright = m_Data->m_Frameslist[ilistLevel][fidx].fEndproj;
				float tcur = m_Data->m_Frameslist[ilistLevel][fidx].fTime;
				//float trange = m_Data->m_Frameslist[ilistLevel][fidx].fprojRange*0.5f;
				//float trange = max(tcur-tleft,1e-3);
				//trange = max(trange, tright- tcur-1);
				cout << tleft << " <-- " << tcur << " -->" << tright << " ||" << trange << endl;
				for (int ik = tleft; ik <= tright; ik++)
				{
					float tdelta = (ik- tcur) / trange;// [-1, 1]
					cout << ik<<"  "<< tdelta <<" ";
					tdelta = min(max(-1.0f, tdelta), 1.0f);
					cout << "tdelta after:" << tdelta << endl;
					WarpAndProjct(ik, fidx, m_voldims, tdelta, aabb, outloop, wpIter, _Base, ilistLevel);
				}
				cout << endl;


			}
		}
		



	} // outloop
	
	

	// init for next level 
	// init index relation [4i]=[2i]=[i]

	if (ilistLevel<nlevels - 1) {
	for (int ifs = 0; ifs < m_nframes; ifs++)
	{
	
		//advection::advect(dimension, aabb, m_Data->m_Frameslist[max(0, iframes - 2)].fFlow, deltaVol, -rtime);
		m_Data->m_Frameslist[ilistLevel + 1][2 * ifs].fVolume = m_Data->m_Frameslist[ilistLevel][ifs].fVolume;
		//m_Data->m_Frameslist[ilistLevel + 1][2 * ifs + 1].fVolume = 
		//advection::advect(m_voldims, aabb, m_Data->m_Frameslist[ilistLevel][ifs].fFlow, m_Data->m_Frameslist[ilistLevel][ifs].fVolume, +0.5f);

		m_Data->m_Frameslist[ilistLevel + 1][2 * ifs + 1].fVolume = m_Data->m_Frameslist[ilistLevel][ifs].fVolume;
		//advection::advect(m_voldims, aabb, m_Data->m_Frameslist[ilistLevel][ifs].fFlow, m_Data->m_Frameslist[ilistLevel][ifs].fVolume, +0.5f);
		// end 


	}


}

//// for  plume 150
//if (ilistLevel == nlevels - 2)
//{
//	cout << "The amount of frames: " << m_nframes << endl;
//
//	for (int ifs = 0; ifs < m_Data->m_fInLevel[ilistLevel+1]; ifs++)
//	{
//		float fi = ifs*1.0f;
//		int indx = min((int)ceilf(fi / 1.25f), m_nframes - 1);
//		//cout << ifs << " :" << indx << endl;;
//		m_Data->m_Frameslist[ilistLevel + 1][ifs].fVolume = m_Data->m_Frameslist[ilistLevel][indx].fVolume;
//		m_Data->m_Frameslist[ilistLevel + 1][ifs].fVolume.save_analyze(str_format("TestFinalLevel/f%03d.hdr", ifs).c_str());
//	}
//	//system("pause");
//
//}

//m_Data->m_Frameslist[ilistLevel + 1][29].fVolume.save_tiff("The29.tiff");
//system("pause");
}// ilistLevel
	exit(0);
	return;
	//Image_date
	/*	delete []Image_data;*/
}


//void PSART::AdaptiveWarp(int iframes, int level, int *dimension, int *projsize, int nb_warp, int outloop, float *aabb)
//{
//	pY_Adaptive.setZero();
//	//PRO_WIDTH = projsize[0];
//	//PRO_HEIGHT = projsize[1];
//	//cout<<"Begin---"<<NO_SART<<"th iteration"<<endl;
//	cout << "Begin---" << idxAdaptive[iframes] << "th iteration for " << iframes << " frame AdaptiveWarp" << endl;
//	image deltaVol = image(dimension[0], dimension[1], dimension[2], 1, 0.0f);
//
//	//	cout<<"============="<<timeofPO++<<"  times in Proximal Operation iteration =============="<<endl;
//	int kkk = 0;
//
//	//////////////////////////////////////////////////////////////////////////
//	// init volume data as U in  SART 
//	//////////////////////////////////////////////////////////////////////////
//
//	// GM changes 5/6/2018 for ContinuousTimeTomo
//	//int dimension[3];
//	////float spacings[3];
//	//m_Data->GetDimensions(dimension);
//
//	int Yindex = dimension[0];
//	int Zindex = dimension[0] * dimension[1];
//
//	//cout<<"==============Begin the  SART algorithm=============="<<endl;
//	//time_t start, end;
//	//time(&start);
//	//int Noiter=ITERATIONS;
//
//	//if(CP_iter!=0)
//	//Noiter=ITERATIONS/2;
//	//ITERATIONS
//	for (int iter = 1;iter <= nb_warp;iter++)
//	{
//		deltaVol.fill(0.0f);
//
//		//cout<<"Checking1 :Factor"<<Factor<<endl;
//		std::srand(unsigned(std::time(0)));
//		std::vector<int> myvector;
//
//		// set some values:  // GM CTT: change to myvector.push_back(i+startproj)
//		//for (int i = 0; i < PRO_NO; ++i) myvector.push_back(i); // 1 2 3 4 5 6 7 8 9
//
//		//// using built-in random generator:
//		//std::random_shuffle(myvector.begin(), myvector.end());
//		//// using myrandom:
//		//std::random_shuffle(myvector.begin(), myvector.end(), myrandom);
//
//		//for (int i = 1;i <= PRO_NO;i++)
//		int _begin = m_Data->m_Frameslist[iframes].fStartproj;
//		int _end = m_Data->m_Frameslist[iframes].fEndproj;
//		for (int i = _begin;i <= _end;i++)
//		{
//			cout << i << "...";
//			float rtime = m_Data->m_ProjIndex[i].rElapsed;
//
//			kkk++;
//			m_NeedCalculateMatrix = true;
//			if (m_NeedCalculateMatrix)
//			{
//				//this->m_Rotation=  pQuaternion[i-1]; 
//
//				//this->m_Rotation = pQuaternion[iframes][myvector[i - 1]];
//				this->m_Rotation = pQuatProjs[i];
//				//cout << "i: " << i << " : " << myvector[i - 1] << "	 " << endl;;
//				/*	cout << "Start to end  projection: " << m_Data->m_Frameslist[iframes].fStartproj << " " << m_Data->m_Frameslist[iframes].fEndproj << endl;
//				cout << " Projection index " << i <<  " From range[ " << m_Data->m_Frameslist[iframes].fStartproj << "," << m_Data->m_Frameslist[iframes].fEndproj <<
//				" ] for key frame: " << iframes << "at time step:" << m_Data->m_Frameslist[iframes].fTime << endl;
//				*/
//
//				//m_RayProducer->m_angle=myvector[i-1]*(Views);
//				_calculateMatrix();
//				m_NeedCalculateMatrix = false;
//			}
//			glMatrixMode(GL_PROJECTION);
//			glLoadMatrixf(m_ProjectionMatrix);
//
//			glMatrixMode(GL_MODELVIEW);
//			glLoadMatrixf(m_ViewMatrix);
//			//	cout<<m_Data->GetProjectionMode()<<endl;
//
//			//	cout << "Check 0" << endl;
//			////////////////////////////////////////////////////////////////////////////////////////////////////		
//			///==================================ForwardProjection==============================================
//			////////////////////////////////////////////////////////////////////////////////////////////////////
//			//RenderContent(iframes, level);
//			//RenderContent(m_Data->m_Frameslist[iframes].fVolume, level);
//			RenderContent(advection::advect(dimension, aabb, m_Data->m_Frameslist[max(0, iframes - 1)].fFlow, m_Data->m_Frameslist[iframes].fVolume, rtime), level);
//
//
//			////////////////////////////////////////////////////////////////////////////////////////////////////	
//			///==================================Corrections==============================================//////
//			////////////////////////////////////////////////////////////////////////////////////////////////////			
//
//			//cout << "Check 1" << endl;
//			//#pragma omp parallel for
//			//for (int k=0;k<PRO_WIDTH*PRO_HEIGHT;k++)
//			//{
//#pragma omp parallel for
//			for (int n = 0;n < PRO_HEIGHT;n++)
//			{
//				for (int m = 0;m < PRO_WIDTH;m++)
//				{
//
//					//int k = m + n*PRO_WIDTH;
//					int k = idx(m, n, PRO_WIDTH);
//					//m_RayProducer->m_Correction[k] = (float)(Factor*Projslist[iframes](m, n, myvector[i - 1]) - Factor * m_RayProducer->GetImageData()[k]) / (float)(Factor*m_RayProducer->GetImageRaylength()[k]);; //weight 为什么如此大？
//
//					//m_RayProducer->m_Correction[k] = (float)(projslevel[level][iframes](m, n, myvector[i - 1]) -  m_RayProducer->GetImageData()[k]) / (float)(m_RayProducer->GetImageRaylength()[k]);; //weight 为什么如此大？
//					//m_RayProducer->m_Correction[k] = (float)(projslevel[level][iframes](m, n, i ) - m_RayProducer->GetImageData()[k]) / (float)(m_RayProducer->GetImageRaylength()[k]);; //weight 为什么如此大？
//					m_RayProducer->m_Correction[k] = (float)(Factor*projslevel[level][i](m, n, 0) - Factor * m_RayProducer->GetImageData()[k] - pY_Adaptive(0, k)) / (float)(Factor*m_RayProducer->GetImageRaylength()[k] + 1.0f);; //weight 为什么如此大？
//					//m_RayProducer->m_Correction[k] = (float)(Factor * projslevel[level][i](m, n, 0) - Factor * m_RayProducer->GetImageData()[k] - pY_warp(0, k)) / (float)(Factor*m_RayProducer->GetImageRaylength()[k] + 1.0f);; //weight 为什么如此大？
//																																																								  //m_RayProducer->m_Correction[k] = (float)(Factor * projslevel[level][i](m, n, 0) - Factor * m_RayProducer->GetImageData()[k] - pY_warp(0, k)) / (float)(Factor*m_RayProducer->GetImageRaylength()[k] + 1.0f);; //weight 为什么如此大？
//
//																																																									//cout << PRO_WIDTH << " " << PRO_HEIGHT << m_RayProducer->m_Correction[k] << endl;
//																																																									//}
//				}
//			}
//			//update pY for PSART algorithm
//			//////////////////////////////////////////////////////////////////////////
//#pragma omp parallel for
//			for (int k = 0;k < PRO_WIDTH*PRO_HEIGHT;k++)
//			{
//				//pY(i-1,k)+=m_alpha*m_RayProducer->m_Correction[k];
//				pY_Adaptive(0, k) += m_alpha * m_RayProducer->m_Correction[k];
//
//				//pY(i-1,k)+=Alpha*m_RayProducer->m_Correction[k];
//				//pY[i-1][k]+=MyAlpha*m_RayProducer->m_Correction[k];
//			}
//
//			////backprojection!
//			//	cout << "======  Voxelbased ==============" << endl;
//#pragma omp parallel for
//			for (int z = m_Data->m_roi[4];z < dimension[2] - m_Data->m_roi[5];z++)
//			{
//				for (int y = m_Data->m_roi[2];y < dimension[1] - m_Data->m_roi[3];y++)
//				{
//					for (int x = m_Data->m_roi[0];x < dimension[0] - m_Data->m_roi[1];x++)
//					{
//
//						/*		if (x<cut[0] || x>dimension[0] - cut[1] || y<cut[2] || y>dimension[1] - cut[3] || z<cut[4] || z>dimension[2] - cut[5])
//						{
//						m_Data->m_Volumelist[iframes](x, y, z) =0.0f;
//						continue;
//						}*/
//						/*		if (!imROI(x, y, z))
//						{
//						m_Data->m_Frameslist[iframes].fVolume(x, y, z) = 0.0f;
//						continue;
//						}*/
//						//int pos=x+y*Yindex+z*Zindex;
//						int pos = idx3D(x, y, z, dimension[0], dimension[1]);
//						Vector4 vec(x, y, z);
//						Vector4 det_pix = m_RayProducer->m_model_view*m_RayProducer->m_VolumeToModelMatrix*vec;
//						float _px = m_halfFrustumWidth - (det_pix[0] / det_pix[2])*m_sdd;
//						float _py = m_halfFrustumHeight - (det_pix[1] / det_pix[2])*m_sdd;
//						_px = _px / m_ds;
//						_py = _py / m_ds;
//						int px = (int)_px;
//						int py = (int)_py;
//
//						//cout << px << " " << py << endl;
//						/*		if (px<boundx1 || px>PRO_WIDTH - boundx2 || py<boundy1 || py>PRO_HEIGHT - boundy2)
//						{
//						m_Data->m_Volumelist[iframes](x, y, z) = 0.0f;
//						continue;
//						}*/
//						float dx = _px - px;
//
//						float dy = _py - py;
//
//						//float rx=1-dx;
//
//						//float ry=1-dy;
//
//						////bicubic interpolation
//						if (px < 2 || px >= PRO_WIDTH - 2 || py < 2 || py >= PRO_HEIGHT - 2)
//							continue;
//						//if(px<2||px>=PRO_WIDTH-2||py<2||py>=PRO_HEIGHT-2)
//						//	continue;
//						////// begin bicubic interpolation
//
//						float *pC = m_RayProducer->m_Correction;
//						pC += py * PRO_WIDTH + px;
//						//double pimg[4][4] = {1,2,3,4,1,10,1000,4,1,100,10000,4,1,2,3,4};
//
//						//if()
//						float pimg[4][4] = { pC[-1 - iminY],pC[-1],pC[-1 + iminY],pC[-1 + 2 * iminY],
//
//							pC[-iminY],pC[0],pC[iminY],pC[2 * iminY],
//
//							pC[1 - iminY],pC[1],pC[1 + iminY],pC[1 + 2 * iminY],
//
//							pC[2 - iminY],pC[2],pC[2 + iminY],pC[2 + 2 * iminY]
//
//						};
//						/*		cout << PRO_WIDTH << " " << PRO_HEIGHT << " " << pC[-1 - iminY] <<" "
//
//						<< pC[-1] <<" "<< pC[-1 + 2 * iminY]<<" "<< pC[2 - iminY]<<" "<< pC[2 + 2 * iminY]<< endl;
//						*/
//
//
//						float gray_value = bicubicInterpolate(pimg, dx, dy);
//						//m_Data->m_Data[pos]+=m_alpha*gray_value;
//						//	Volumelist[iframes](x, y, z) += 0.5*gray_value;
//
//						deltaVol(x, y, z) += gray_value;
//						//	m_Data->m_Frameslist[iframes].fVolume(x, y, z) = max(min(G_Max, m_Data->m_Frameslist[iframes].fVolume(x, y, z)), 0);
//						//m_Data->m_VolumeLevel[level][iframes](x, y, z) += m_alpha * gray_value;
//						//m_Data->m_VolumeLevel[level][iframes](x, y, z) = max(min(G_Max, m_Data->m_VolumeLevel[level][iframes](x, y, z)), 0);
//
//					}
//
//				}
//
//				//generate a ray from ray source to voxel center shooting the detector
//
//			}
//
//
//			m_Data->m_Frameslist[iframes].fVolume += m_alpha * advection::advect(dimension, aabb, m_Data->m_Frameslist[max(0, iframes - 2)].fFlow, deltaVol, -rtime);
//			//m_Data->m_Frameslist[iframes].fVolume = max(m_Data->m_Frameslist[iframes].fVolume, 0.0f);
//
//			m_Data->m_Frameslist[iframes].fVolume.cut(0.0, G_Max);
//
//			//cout << "Check 3" << endl;
//			m_RayProducer->Initial();
//			m_Data->SetProjectionMode(0);
//		}
//
//		cout << endl;
//	}
//
//
//
//	//cout << "Check 4" << endl;
//	int dimenst[3];
//	m_Data->GetDimensions(dimenst);
//	char Mybuf[10];
//	//sprintf(Mybuf, "%d", NO_SART);
//	sprintf(Mybuf, "%d", idxAdaptive[iframes]);
//	char Mybuf3[10];
//	sprintf(Mybuf3, "%03d", iframes);
//	char Mybuf4[10];
//	sprintf(Mybuf4, "%03d", level);
//	char Mybuf5[10];
//	sprintf(Mybuf5, "%03d", outloop);
//
//
//	string Myb = Mybuf;
//	string Myb3 = Mybuf3;
//	string Myb4 = Mybuf4;
//	string Myb5 = Mybuf5;
//	string _output = pOutput;
//	//string Mywritefile=output+Myb+outputformat;
//	//int pos=_output.find(".mha");
//	//if(pos>1000) pos=10;
//	//string Mywritefile=_output.substr(0,pos)+"_"+Myb+"iter.mha";
//
//	//string Mywritefile=_output+"-"+Myb3+"-frame-"+Myb+".mha";
//	string Mywritefile2 = _output + "-" + Myb3 + "-frame-" + Myb4 + "-level-" + Myb + "iter_VolOpt-" + Myb5 + "_outloop.tiff";
//	cout << Mywritefile2 << endl;
//	//char* MyOutputfile = strdup(Mywritefile.c_str());
//	char* MyOutputfile2 = strdup(Mywritefile2.c_str());
//
//	int dims2[3];
//	float spacs2[3];
//	////if((CP_iter%4==0))
//	////{
//
//	//	FILE *fp5;
//	//	fp5=fopen(MyOutputfile,"wb");
//
//	//	m_Data->GetDimensions(dims2);
//	//	m_Data->GetSpacings(spacs2);
//	//	fwrite(pOutputFile,sizeof(char),len,fp5);
//	//	//fwrite(m_Data->GetData(),sizeof(float),dims2[0]*dims2[1]*dims2[2],fp5);
//
//	//	fwrite(Volumelist[iframes], sizeof(float), dims2[0] * dims2[1] * dims2[2], fp5);
//	//	fclose(fp5);
//	//
//	//	//}
//
//#pragma omp parallel for
//	for (int z = 0;z < dimension[2] - 0;z++)
//	{
//		for (int y = 0;y < dimension[1] - 0;y++)
//		{
//			for (int x = 0;x < dimension[0] - 0;x++)
//			{
//
//				/*	if (x<cut[0] || x>dimension[0] - cut[1] || y<cut[2] || y>dimension[1] - cut[3] || z<cut[4] || z>dimension[2] - cut[5])
//				{
//				m_Data->m_Volumelist[iframes](x, y, z) = 0.0f;
//				continue;
//				}*/
//				//m_Data->m_Volumelist[iframes](x, y, z) *= imROI(x, y, z);
//				m_Data->m_Frameslist[iframes].fVolume(x, y, z) = max(m_Data->m_Frameslist[iframes].fVolume(x, y, z)*imROI(x, y, z), 0.0f);
//
//			}
//		}
//	}
//
//	//if(level==0)
//	//	m_Data->m_Frameslist[iframes].fVolume = m_Data->m_Frameslist[iframes].fVolume.get_normalize(0.0, 50.0);
//	//if ((idxSART[iframes] % 2) == 1)
//	//m_Data->m_VolumeLevel[level][iframes].save_tiff(MyOutputfile2);
//	//if (iframes == 2)
//		m_Data->m_Frameslist[iframes].fVolume.save_tiff(MyOutputfile2);
//
//	//m_Data->m_Frameslist[iframes].fVolume.normalize(0.0, G_Max).save_tiff(MyOutputfile2);;
//
//	//m_Data->m_Volumelist[iframes].save_tiff(MyOutputfile2);
//
//	//m_Data->m_Volumelist[iframes].save_analyze(MyOutputfile2);
//	//Volumelist[iframes].save_tiff(MyOutputfile2, 0);
//
//	//NO_SART++;
//	cout << "end---" << idxAdaptive[iframes] << "th iteration for AdaptiveWarp" << iframes << " frame" << endl;
//	cout << "Finish " << iframes << " frame " << idxAdaptive[iframes] << "times primal part calculation." << endl;
//	idxAdaptive[iframes]++;
//
//	//cout<<endl;
//
//}
//

void PSART::WarpAndProjct(int iproj, int iframes, int *dimension, float rtime, float *aabb, int outloop, int warpiter, char * _file_header, int flevel)
{
	pY_warp.setZero();
	//advection::advect(dim, aabb, m_Data->m_Frameslist[rFrameidx].fFlow,m_Data->m_Frameslist[rFrameidx].fVolume, rTime)
	int kkk = 0;
	image deltaVol = image(dimension[0], dimension[1], dimension[2], 1, 0.0f);
	int Yindex = dimension[0];
	int Zindex = dimension[0] * dimension[1];
	//cout << "Test 1" << endl;
	for (int iter = 1;iter < 2;iter++)
	{
		deltaVol = 0.0;
		//cout<<"Checking1 :Factor"<<Factor<<endl;
		//std::srand(unsigned(std::time(0)));
		std::vector<int> myvector;


		int proId = iproj;
		//RotateVolume(i);
		RotateVolume(proId);

		kkk++;
		if (m_NeedCalculateCamera)
		{
			_calculateCamera();
			m_NeedCalculateCamera = false;
		}
		if (m_NeedCalculateMatrix)
		{
			_calculateMatrix();
			m_NeedCalculateMatrix = false;
		}
		glMatrixMode(GL_PROJECTION);
		glLoadMatrixf(m_ProjectionMatrix);

		glMatrixMode(GL_MODELVIEW);
		glLoadMatrixf(m_ViewMatrix);



		
		//	cout<<m_Data->GetProjectionMode()<<endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////		
		///==================================ForwardProjection==============================================
		//advection::advect(dimension, aabb, m_Data->m_Frameslist[max(0, iframes - 1)].fFlow, m_Data->m_Frameslist[iframes].fVolume, rtime)
		RenderContent(iframes,true, iproj, rtime,flevel);

		///==================================Corrections==============================================//////
		////////////////////////////////////////////////////////////////////////////////////////////////////			

		//cout << "Check 1" << endl;
		//#pragma omp parallel for
		//for (int k=0;k<PRO_WIDTH*PRO_HEIGHT;k++)
		//{
#pragma omp parallel for
		for (int n = 0;n < PRO_HEIGHT;n++)
		{
			for (int m = 0;m < PRO_WIDTH;m++)
			{

				//int k = m + n*PRO_WIDTH;
				int k = idx(m, n, PRO_WIDTH);
				
				m_RayProducer->m_Correction[k] = (float)(Factor * projslevel[0][proId](m, n, 0) - Factor * m_RayProducer->GetImageData()[k] - pY_warp(0, k)) / (float)(Factor*m_RayProducer->GetImageRaylength()[k] + 1.0f);; //weight 为什么如此大？
			
			}
		}
		//update pY for PSART algorithm
		//////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
		for (int k = 0;k < PRO_WIDTH*PRO_HEIGHT;k++)
		{
			//pY(i-1,k)+=m_alpha*m_RayProducer->m_Correction[k];
			pY_warp(0, k) += m_alpha * m_RayProducer->m_Correction[k];

			//pY(i-1,k)+=Alpha*m_RayProducer->m_Correction[k];
			//pY[i-1][k]+=MyAlpha*m_RayProducer->m_Correction[k];
		}


		////backprojection!

		//cout << "Check 2" << endl;
		//if (Backprojector == Voxelbased)
		//{

			//	cout << "======  Voxelbased ==============" << endl;
#pragma omp parallel for
		for (int z = vh_z1;z < dimension[2] - vh_z2;z++)
		{
			for (int y = vh_y1;y < dimension[1] - vh_y2;y++)
			{
				for (int x = vh_x1;x < dimension[0] - vh_x2;x++)
				{

					/*		if (x<cut[0] || x>dimension[0] - cut[1] || y<cut[2] || y>dimension[1] - cut[3] || z<cut[4] || z>dimension[2] - cut[5])
					{
					m_Data->m_Volumelist[iframes](x, y, z) =0.0f;
					continue;
					}*/
					//int pos=x+y*Yindex+z*Zindex;

					int pos = idx3D(x, y, z, dimension[0], dimension[1]);
					Vector4 vec(x, y, z);
					Vector4 det_pix = m_RayProducer->m_model_view*m_RayProducer->m_VolumeToModelMatrix*vec;
					float _px = m_halfFrustumWidth - (det_pix[0] / det_pix[2])*m_sdd;
					float _py = m_halfFrustumHeight - (det_pix[1] / det_pix[2])*m_sdd;
					_px = _px / m_ds;
					_py = _py / m_ds;
					int px = (int)_px;
					int py = (int)_py;

					//cout << px << " " << py << endl;
					/*		if (px<boundx1 || px>PRO_WIDTH - boundx2 || py<boundy1 || py>PRO_HEIGHT - boundy2)
					{
					m_Data->m_Volumelist[iframes](x, y, z) = 0.0f;
					continue;
					}*/
					float dx = _px - px;

					float dy = _py - py;

					//float rx=1-dx;

					//float ry=1-dy;

					////bicubic interpolation
					if (px < 2 || px >= PRO_WIDTH - 2 || py < 2 || py >= PRO_HEIGHT - 2)
						continue;
					//if(px<2||px>=PRO_WIDTH-2||py<2||py>=PRO_HEIGHT-2)
					//	continue;
					////// begin bicubic interpolation

					float *pC = m_RayProducer->m_Correction;
					pC += py * PRO_WIDTH + px;
					//double pimg[4][4] = {1,2,3,4,1,10,1000,4,1,100,10000,4,1,2,3,4};

					//if()
					float pimg[4][4] = { pC[-1 - iminY],pC[-1],pC[-1 + iminY],pC[-1 + 2 * iminY],

						pC[-iminY],pC[0],pC[iminY],pC[2 * iminY],

						pC[1 - iminY],pC[1],pC[1 + iminY],pC[1 + 2 * iminY],

						pC[2 - iminY],pC[2],pC[2 + iminY],pC[2 + 2 * iminY]

					};
					/*		cout << PRO_WIDTH << " " << PRO_HEIGHT << " " << pC[-1 - iminY] <<" "

					<< pC[-1] <<" "<< pC[-1 + 2 * iminY]<<" "<< pC[2 - iminY]<<" "<< pC[2 + 2 * iminY]<< endl;
					*/


					float gray_value = bicubicInterpolate(pimg, dx, dy);
					//m_Data->m_Data[pos]+=m_alpha*gray_value;
					//	Volumelist[iframes](x, y, z) += 0.5*gray_value;

		
					deltaVol(x, y, z) += gray_value;
					

				}

			}

			//generate a ray from ray source to voxel center shooting the detector

		}


		//m_Data->m_Frameslist[iframes].fVolume += m_alpha * advection::advect(dimension, aabb, m_Data->m_Frameslist[max(0, iframes - 2)].fFlow, deltaVol, -rtime);
		m_Data->m_Frameslist[flevel][iframes].fVolume += m_alpha * advection::advect(dimension, aabb, m_Data->m_Frameslist[flevel][max(0, iframes - 1)].fFlow, deltaVol, -rtime);
		//m_Data->m_Frameslist[iframes].fVolume = max(m_Data->m_Frameslist[iframes].fVolume, 0.0f);

		m_Data->m_Frameslist[flevel][iframes].fVolume.cut(0.0, G_Max);

		deltaVol.assign();

		//cout << "Check 3" << endl;
		m_RayProducer->Initial();
		m_Data->SetProjectionMode(0);
		//}

	//	cout << endl;
	}


	

	int dims2[3];
	float spacs2[3];
#pragma omp parallel for
	for (int z = 0;z < dimension[2] - 0;z++)
	{
		for (int y = 0;y < dimension[1] - 0;y++)
		{
			for (int x = 0;x < dimension[0] - 0;x++)
			{

		
				m_Data->m_Frameslist[flevel][iframes].fVolume(x, y, z) = max(m_Data->m_Frameslist[flevel][iframes].fVolume(x, y, z)*imROI(x, y, z), 0.0f);


			}
		}
	}



	//cout << "Test 3" << endl;
	if(warpiter==0)
		m_Data->m_Frameslist[flevel][iframes].fVolume.save_analyze(str_format("%s/f%03d.p%03d.warp%03d.out%02d.flevel%02d.hdr", _file_header,  iframes, iproj,  warpiter,outloop,flevel).c_str());


}


void PSART::SART_Proximal_Operator(int iframes,  int nb_sart,int outloop, int stIter, int flevel )
{
	pY.setZero();

	//PRO_WIDTH = projsize[0];
	//PRO_HEIGHT = projsize[1];
	//cout<<"Begin---"<<NO_SART<<"th iteration"<<endl;
	cout << "Begin---" << idxSART_Proximal[iframes] << "th iteration for " << iframes << " frame" << endl;

	//	cout<<"============="<<timeofPO++<<"  times in Proximal Operation iteration =============="<<endl;
	int kkk = 0;

	//////////////////////////////////////////////////////////////////////////
	// init volume data as U in  SART 
	//////////////////////////////////////////////////////////////////////////

	for (int iter = 1;iter <= nb_sart;iter++)
	{
		//cout<<"Checking1 :Factor"<<Factor<<endl;
		std::srand(unsigned(std::time(0)));
		std::vector<int> myvector;

		// set some values:  

	//	int _begin = m_Data->m_Frameslist[iframes].fStartproj+1;
		int _begin = m_Data->m_Frameslist[flevel][iframes].fStartproj ;
		int _end = m_Data->m_Frameslist[flevel][iframes].fEndproj;

		int first = _begin;
		int base = first;
		int end = _end;
		int ind = 0;

		int num = end - first + 1;

		int *front = new int[num];
		int *back = new int[num];
		int *result = new int[num];

		for (int i = first; i < end + 1; i++)
		{
			front[ind] = i;
			back[ind] = end - ind;
			ind++;
		}

		for (int i = 0; i < num; i++)
		{
			if (2 * i<num)
				result[2 * i] = front[i];
			else
				break;
		}

		for (int i = 0; i < num; i++)
		{
			if ((2 * i + 1)<num)
				result[2 * i + 1] = back[i];
			else
				break;
		}




		//for (int i = _begin-1;i <= _end;i++)
		for (int i = _begin ;i <= _end;i++)
		{
			cout << i << "...";// 
			
			

			RotateVolume(result[i- _begin ]);

			kkk++;
			if (m_NeedCalculateCamera)
			{
				_calculateCamera();
				m_NeedCalculateCamera = false;
			}
			if (m_NeedCalculateMatrix)
			{
				_calculateMatrix();
				m_NeedCalculateMatrix = false;
			}
			glMatrixMode(GL_PROJECTION);
			glLoadMatrixf(m_ProjectionMatrix);

			glMatrixMode(GL_MODELVIEW);
			glLoadMatrixf(m_ViewMatrix);
			//	cout<<m_Data->GetProjectionMode()<<endl;

			
			////////////////////////////////////////////////////////////////////////////////////////////////////		
			///==================================ForwardProjection==============================================
			////////////////////////////////////////////////////////////////////////////////////////////////////
			//RenderContent(iframes, level);
			RenderContent(iframes,0,0,0,flevel);


			////////////////////////////////////////////////////////////////////////////////////////////////////	
			///==================================Corrections==============================================//////
			////////////////////////////////////////////////////////////////////////////////////////////////////			

		

#pragma omp parallel for
			for (int n = 0;n < PRO_HEIGHT;n++)
			{
				for (int m = 0;m < PRO_WIDTH;m++)
				{

					//int k = m + n*PRO_WIDTH;
					int k = idx(m, n, PRO_WIDTH);

					m_RayProducer->m_Correction[k] = (float)(Factor*projslevel[0][result[i - _begin]](m, n, 0) - Factor * m_RayProducer->GetImageData()[k] - pY(result[i - _begin], k)) / (float)(Factor*m_RayProducer->GetImageRaylength()[k] + 1.0f);; 

																																			
				}
			}



////			//update pY for PSART algorithm
////			//////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
			for (int k = 0;k < PRO_WIDTH*PRO_HEIGHT;k++)
			{

				pY(result[i - _begin], k) += m_alpha * m_RayProducer->m_Correction[k];
			
			}

#pragma omp parallel for
			for (int z = m_Data->m_roi[4];z < m_voldims[2] - m_Data->m_roi[5];z++)
			{
				for (int y = m_Data->m_roi[2];y < m_voldims[1] - m_Data->m_roi[3];y++)
				{
					for (int x = m_Data->m_roi[0];x < m_voldims[0] - m_Data->m_roi[1];x++)
					{

					
								//int pos=x+y*Yindex+z*Zindex;
						int pos = idx3D(x, y, z, m_voldims[0], m_voldims[1]);
						Vector4 vec(x, y, z);
						Vector4 det_pix = m_RayProducer->m_model_view*m_RayProducer->m_VolumeToModelMatrix*vec;
						float _px = m_halfFrustumWidth - (det_pix[0] / det_pix[2])*m_sdd;
						float _py = m_halfFrustumHeight - (det_pix[1] / det_pix[2])*m_sdd;
						_px = _px / m_ds;
						_py = _py / m_ds;
						int px = (int)_px;
						int py = (int)_py;

						//cout << px << " " << py << endl;
						/*		if (px<boundx1 || px>PRO_WIDTH - boundx2 || py<boundy1 || py>PRO_HEIGHT - boundy2)
						{
						m_Data->m_Volumelist[iframes](x, y, z) = 0.0f;
						continue;
						}*/
						float dx = _px - px;

						float dy = _py - py;

					
						////bicubic interpolation
						if (px < 2 || px >= PRO_WIDTH - 2 || py < 2 || py >= PRO_HEIGHT - 2)
							continue;
						
						////// begin bicubic interpolation

						float *pC = m_RayProducer->m_Correction;
						pC += py * PRO_WIDTH + px;
						//double pimg[4][4] = {1,2,3,4,1,10,1000,4,1,100,10000,4,1,2,3,4};

						//if()
						float pimg[4][4] = { pC[-1 - iminY],pC[-1],pC[-1 + iminY],pC[-1 + 2 * iminY],

							pC[-iminY],pC[0],pC[iminY],pC[2 * iminY],

							pC[1 - iminY],pC[1],pC[1 + iminY],pC[1 + 2 * iminY],

							pC[2 - iminY],pC[2],pC[2 + iminY],pC[2 + 2 * iminY]

						};
						/*		cout << PRO_WIDTH << " " << PRO_HEIGHT << " " << pC[-1 - iminY] <<" "

						<< pC[-1] <<" "<< pC[-1 + 2 * iminY]<<" "<< pC[2 - iminY]<<" "<< pC[2 + 2 * iminY]<< endl;
						*/


						float gray_value = bicubicInterpolate(pimg, dx, dy);
						
						m_Data->m_Frameslist[flevel][iframes].fVolume(x, y, z) += m_alpha * gray_value;

						m_Data->m_Frameslist[flevel][iframes].fVolume(x, y, z) = max(min(G_Max, m_Data->m_Frameslist[flevel][iframes].fVolume(x, y, z)), 0);

					

					}

				}

				//generate a ray from ray source to voxel center shooting the detector

			}


			//cout << "Check 3" << endl;
			m_RayProducer->Initial();
			m_Data->SetProjectionMode(0);
		}

		cout << endl;
	}




	int dims2[3];
	float spacs2[3];


#pragma omp parallel for
	for (int z = 0;z < m_voldims[2] - 0;z++)
	{
		for (int y = 0;y < m_voldims[1] - 0;y++)
		{
			for (int x = 0;x < m_voldims[0] - 0;x++)
			{

					//m_Data->m_Volumelist[iframes](x, y, z) *= imROI(x, y, z);
				m_Data->m_Frameslist[flevel][iframes].fVolume(x, y, z) = max(m_Data->m_Frameslist[flevel][iframes].fVolume(x, y, z)*imROI(x, y, z), 0.0f);

			}
		}
	}

	//cout << "OK new" << endl;
	char out_dir[1024];
	sprintf(out_dir, pOutput);

	
	int ft = m_Data->m_Frameslist[flevel][iframes].fTime;
	int range = m_Data->m_Frameslist[flevel][iframes].fprojRange;
	//if (!(stIter%3))
		m_Data->m_Frameslist[flevel][iframes].fVolume.save_analyze(str_format("%s.f%03d.out%02d.st%02d.flevel%02d.hdr", out_dir, iframes, outloop,stIter, flevel).c_str());


	cout << "end---" << idxSART_Proximal[iframes] << "th iteration for SART Proximal Operators" << iframes << " frame" << endl;
	cout << "Finish " << iframes << " frame " << idxSART_Proximal[iframes] << "times primal part calculation." << endl;
	idxSART_Proximal[iframes]++;

	//cout<<endl;

}





//settings for OPENGL
//void PSART::RenderContent(int frame, int level)
void PSART::RenderContent(int frame, bool isWarpped, int ProIdx, float rtime, int flevel)

{


	glClearColor(m_BackgroundColor[0], m_BackgroundColor[1], m_BackgroundColor[2], 1.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	m_RayProducer->SetViewMatrix(this->GetViewMatrix());
	m_RayProducer->SetProjectionMatrix(this->GetProjectionMatrix());


	if (_renderCoreFuc(frame, isWarpped, ProIdx, rtime, flevel)) m_RayProducer->RenderImage(); //attention!

}


void PSART::GetSize3D(float size[3])
{
	int dims[3];
	float spacings[3];
	if (m_Data)
	{
		m_Data->GetDimensions(dims);
		m_Data->GetSpacings(spacings);
	}

	size[0] = dims[0] * spacings[0];
	size[1] = dims[1] * spacings[1];
	size[2] = dims[2] * spacings[2];
}

#include "SingleScalarVolume.core.h"
#include "RayCaster.core.h"
#include "TypicalRayProducer.core.h"


//#include "GridPartition.core.h"
#include <iostream>
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string>
using namespace std;





// static float * FinalVolumeData= new float[PRO_WIDTH*PRO_HEIGHT] ;
static VolumeData *vd;//

//bool PSART::_renderCoreFuc(int frame, int level)
bool PSART::_renderCoreFuc(int frame, bool isWarpped , int ProIdx, float rtime, int flevel)
{


	////m_Data->SetData(VoxelValue);
	//if (!SingleScalarVolumePrepare(m_Data, m_RayProducer, frame, level)) return false;
	if (!SingleScalarVolumePrepare(m_Data, m_RayProducer, frame, isWarpped, ProIdx,  rtime, flevel)) return false;


	if (!RayCasterPrepare(m_RayCaster, m_voldims)) return false;


	GenerateImage_GetCorrections(m_RayProducer);


	return true;





}
//
bool PSART::RotateVolume(int i)
{
	this->m_NeedCalculateMatrix = true;

	this->m_NeedCalculateCamera = true;

	//this->_makeReference();

	double _Index = m_RadPerPixel;


	double angle = inputdegree[i] * 3.1415926 / 180;

	Vector4 v;

	v[0] = 0;  

	v[1] = 1;

	v[2] = v[3] = 0.0f;

	if (!IS_COUNTERCLOCK)

		v[1] *= -1.0f;
	double length = v.Length();
	v *= 1.0f / length;
	//length*
	this->m_Rotation = Quaternion(v[0], v[1], v[2], length*angle)*this->m_OldRotation;


	index_of_project++;
	idx_proj++;
	return true;
}

//

void PSART::InitViews()
{
	// important, need to change flexible
	//int scannerRound = 6;
	// try again
	for (int i = 0; i < scannerRound; i++)
	{
		index_of_project = 0;
		for (size_t j = 0; j < m_nrounds; j++)
		{


			this->m_Rotation.Identity();
			m_OldRotation = m_Rotation;
			this->m_NeedCalculateMatrix = true;
			this->m_NeedCalculateCamera = true;
			this->_makeReference();
			double angle = this->m_startDegree[i][j] * 3.1415926 / 180;
			cout << " angle: " << this->m_startDegree[i][j] << " ";
			//double angle = this->m_RoundstartDegree[j] * 3.1415926 / 180;
			//cout << "angle:  "<<angle ;
			Vector4 v;


			v[0] = 0;  //-15 大概旋转了21次每面  -315半面  27 44次
			v[1] = 1;
			v[2] = v[3] = 0.0f;
			if (!IS_COUNTERCLOCK)
				v[1] *= -1.0f;

			double length = v.Length();
			v *= 1.0f / length;
			//length*
			this->m_Rotation = Quaternion(v[0], v[1], v[2], length*angle)*this->m_OldRotation;

			//cout<<"Rotation:"<<v[0]<<"  "<<v[1]<<" "<<v[2]<<" "<<" Angle "<<length*angle<<endl;
			//cout << " " << " Angle " << length*angle << endl;
			//pQuaternion[i][index_of_project] = this->m_Rotation;
			pQuatProjs[idx_proj] = this->m_Rotation;
			//cout << " " << pQuaternion[i][index_of_project];
			index_of_project++;
			idx_proj++;
			//// try again
			//pQuaternion[index_of_project] = this->m_Rotation;
			//index_of_project++;
		//	for (int index = 0;index<PRO_NO - 1;index++)
			for (int index = 0;index < projsEachRound - 1;index++)
			{
				//RotateVolume();

				if (m_NeedCalculateCamera)
				{
					_calculateCamera();
					m_NeedCalculateCamera = false;
					//	firstCal=false;
					//	cout<<"sdfsdfas"<<endl;

				}
				m_NeedCalculateMatrix = true;
				if (m_NeedCalculateMatrix)
				{
					_calculateMatrix();
					m_NeedCalculateMatrix = false;
				}


				glMatrixMode(GL_PROJECTION);
				glLoadMatrixf(m_ProjectionMatrix);

				glMatrixMode(GL_MODELVIEW);
				glLoadMatrixf(m_ViewMatrix);
				//cout<<endl;
				RotateVolume(i);
			}
			cout << endl;
			cout << endl;
		}

		cout << "NumOf projection test:" << index_of_project << endl;
	
	}

	cout << "Testing for the # of idx_proj  and m_projNo -1:	" << idx_proj << " ---VS--- " << m_projNo << endl;
}




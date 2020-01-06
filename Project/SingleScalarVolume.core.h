// SingleScalarVolume.core.h

#include "VolumeData.h"
#include "TypicalRayProducer.h"
#include "Vector4.h"
#include"advection.h"
#include<iostream>
using namespace std;
CImg<float> idata;
static float *pWeightForVoxels;
static float *pContributions;
static float *pCorrection;
static int m_dims[3];
static float m_spacings[3];
static int m_incs[3];
static bool m_IsBackProjectionMode;

class VolumeAccess
{
public:
	VolumeAccess(){}
	virtual ~VolumeAccess(){}

	// visit sample  points
	virtual float GetScalarValue(int x,int y,int z,float &SumWeight_Ray, int indexOfimage)=0;
	// get the interpolation position
	//virtual float GetScalarValue(Vector4 pos)=0;
	virtual float GetScalarValue(Vector4 pos,float &SumWeight_Ray,int indexOfimage)=0;
	// visit interpolation  kaiser bessel
	virtual float GetScalarValue_Kaiser(Vector4 pos, float theta)=0;
	// get normal
	virtual Vector4 GetNormal(Vector4 pos)=0;
};


//template <typename T>
class T_VolumeAccess : public VolumeAccess
{
private:
	//T* m_data;
	//float *m_data;
	//Cimg<float> idata;

public: 
	//T_VolumeAccess(void* data, const int* dims, const float* spacings, bool IsBackProMode, float * top, float*down, float *correction)

	T_VolumeAccess(CImg<float> data, const int* dims, const float* spacings, bool IsBackProMode, float * top, float*down, float *correction)
	{
		//m_data=(float*)data;
		idata = data;
		pContributions=(float*)top;
		pWeightForVoxels=(float*)down;
		pCorrection=(float*)correction;
		m_IsBackProjectionMode=IsBackProMode;
		m_dims[0]=dims[0];m_dims[1]=dims[1];m_dims[2]=dims[2];
		m_spacings[0]=spacings[0];m_spacings[1]=spacings[1];m_spacings[2]=spacings[2];

		m_incs[0]=1;
		m_incs[1]=m_dims[0];
		m_incs[2]=m_dims[0]*m_dims[1];
	

	}
	virtual ~T_VolumeAccess(){}

	virtual float GetScalarValue(int x,int y,int z,float &SumWeight_Ray, int indexOfimage)
	//virtual float GetScalarValue(int x,int y,int z)
	{
		//cout<<"GetScalarValue is x y z "<<endl;
		//SumWeight_Ray+=1.0f; //weight1

		x=max(min(x,m_dims[0]-1),0);
		y=max(min(y,m_dims[1]-1),0);
		z=max(min(z,m_dims[2]-1),0);

	/*	if(m_IsBackProjectionMode)
		{
				pContributions[x+(y+z*m_dims[1])*m_dims[0]]+=pCorrection[indexOfimage];
				pWeightForVoxels[x+(y+z*m_dims[1])*m_dims[0]]+=1.0f;
		}*/
		//return (float)m_data[x+(y+z*m_dims[1])*m_dims[0]];
		return idata(x,y,z);
	}

	virtual void UpdateContributionSum(int x,int y,int z,int iImag,float w)
	{
	/*	x=max(min(x,m_dims[0]-1),0);
		y=max(min(y,m_dims[1]-1),0);
		z=max(min(z,m_dims[2]-1),0);*/

		pContributions[x+(y+z*m_dims[1])*m_dims[0]]+=pCorrection[iImag]*w;
		pWeightForVoxels[x+(y+z*m_dims[1])*m_dims[0]]+=w;



	}

	//trilinear interpolation
	virtual float GetScalarValue(Vector4 pos,float &SumWeight_Ray, int indexOfimage)
	//virtual float GetScalarValue(Vector4 pos)
	{
		pos[0]=max(min(pos[0],m_dims[0]-1),0);
		pos[1]=max(min(pos[1],m_dims[1]-1),0);
		pos[2]=max(min(pos[2],m_dims[2]-1),0);
		//pos[0] = max(min(pos[0], m_dims[0] - 22), 22);
		//pos[1] = max(min(pos[1], m_dims[1] - 22), 22);
		//pos[2] = max(min(pos[2], m_dims[2] - 22), 22);
		
		float a,b,c,d;
		int ipos[3];
		float fpos1,fpos2;

		ipos[0]=(int)pos[0]; ipos[1]=(int)pos[1]; ipos[2]=(int)pos[2];

			float t0 = pos[0] - ipos[0];
			float t1 = pos[1] - ipos[1];
			float t2 = pos[2] - ipos[2];
			float M_t0 = 1.0f - t0;
			float M_t1 = 1.0f - t1;
			float M_t2 = 1.0f - t2;
			float w[8] = { M_t0*M_t1*M_t2, t0*M_t1*M_t2, M_t0*t1*M_t2, t0*t1*M_t2, M_t0*M_t1*t2, t0*M_t1*t2, M_t0*t1*t2, t0*t1*t2 };



		fpos2=pos[0]-ipos[0];
		fpos1=1.0f-fpos2;
		if (!m_IsBackProjectionMode)
			SumWeight_Ray+=1.0f;
		a= GetScalarValue(ipos[0],ipos[1],ipos[2],SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1],ipos[2],SumWeight_Ray,indexOfimage)*fpos2;
		b= GetScalarValue(ipos[0],ipos[1]+1,ipos[2],SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1]+1,ipos[2],SumWeight_Ray,indexOfimage)*fpos2;
		c= GetScalarValue(ipos[0],ipos[1],ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1],ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos2;
		d= GetScalarValue(ipos[0],ipos[1]+1,ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos1+ GetScalarValue(ipos[0]+1,ipos[1]+1,ipos[2]+1,SumWeight_Ray,indexOfimage)*fpos2;

		fpos2=pos[1]-ipos[1];
		fpos1=1.0f-fpos2;

		a=a*fpos1+b*fpos2;
		c=c*fpos1+d*fpos2;

		fpos2=pos[2]-ipos[2];
		fpos1=1.0f-fpos2;
//
//
//
		//float gw=ptr[0]*M_t0*M_t1*M_t2+ptr[1]*t0*M_t1*M_t2+ptr[m_incs[1]]*M_t0*t1*M_t2+ptr[m_incs[1]+1]*t0*t1*M_t2+
		//	ptr[m_incs[2]]*M_t0*M_t1*t2+ptr[m_incs[2]+1]*t0*t2*M_t1+ptr[m_incs[1]+m_incs[2]]*M_t0*t1*t2+ptr[m_incs[1]+m_incs[2]+1]*t0*t1*t2;
			//int m_loc = ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2];
			//ptr += m_loc;
			if (m_IsBackProjectionMode)
			{
				//if (pos[0] == 2 || (pos[0] == m_dims[0] - 2) ||
				//	pos[1] == 2 || (pos[1] == m_dims[1] - 2) ||
				//	pos[2] == 2 || (pos[2] == m_dims[2] - 2))
				//{
				//	w[0] = w[1] = w[2] = w[3] = w[4] = w[5] = w[6] = w[7] = 0.0f;
				//}

				UpdateContributionSum(ipos[0],ipos[1],ipos[2],indexOfimage,w[0]);
				UpdateContributionSum(ipos[0]+1,ipos[1],ipos[2],indexOfimage,w[1]);
				UpdateContributionSum(ipos[0],ipos[1]+1,ipos[2],indexOfimage,w[2]);
				UpdateContributionSum(ipos[0]+1,ipos[1]+1,ipos[2],indexOfimage,w[3]);
				UpdateContributionSum(ipos[0],ipos[1],ipos[2]+1,indexOfimage,w[4]);
				UpdateContributionSum(ipos[0]+1,ipos[1],ipos[2]+1,indexOfimage,w[5]);
				UpdateContributionSum(ipos[0],ipos[1]+1,ipos[2]+1,indexOfimage,w[6]);
				UpdateContributionSum(ipos[0]+1,ipos[1]+1,ipos[2]+1,indexOfimage,w[7]);



				//float* pCtrb = pContributions;
				////cout<<m_incs[0]<<"		"<<m_incs[1]<<"		"<<"	"<<m_incs[2]<<endl;
				////if(pCorrection[indexOfimage]>2.0)
				////cout<<"pCorrection: "<<  pCorrection[indexOfimage]<<endl;
				////pCtrb += m_loc;
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+0] += (float)(pCorrection[indexOfimage] * w[0]);
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+1] += (float)(pCorrection[indexOfimage] * w[1]);
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1]] += (float)(pCorrection[indexOfimage] * w[2]);
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1] + 1] += (float)(pCorrection[indexOfimage] * w[3]);
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[2]] += (float)(pCorrection[indexOfimage] * w[4]);
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[2] + 1] += (float)(pCorrection[indexOfimage] * w[5]);
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1] + m_incs[2]] += (float)(pCorrection[indexOfimage] * w[6]);
				//pCtrb[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1] + m_incs[2] + 1] += (float)(pCorrection[indexOfimage] * w[7]);
				//if (M_t0*M_t1*M_t2>1.0f || t0*M_t1*M_t2>1.0 || M_t0*t1*M_t2>1.0 || t0*t1*M_t2>1.0 || M_t0*M_t1*t2>1.0 || t0*t2*M_t1>1.0 || M_t0*t1*t2>1.0 || t0*t1*t2>1.0)
				//	cout << M_t0*M_t1*M_t2 << " " << t0*M_t1*M_t2 << "	" << M_t0*t1*M_t2 << "" << t0*t1*M_t2 << " " << M_t0*M_t1*t2 << " " << t0*t2*M_t1 << " " << M_t0*t1*t2 << " " << t0*t1*t2 << endl;
				//
				//float *pweight = pWeightForVoxels;
				////pweight += m_loc;

				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+0] += w[0];
				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+1] += w[1];
				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1]] += w[2];
				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1] + 1] += w[3];
				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[2]] += w[4];
				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[2] + 1] += w[5];
				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1] + m_incs[2]] += w[6];
				//pweight[ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2]+m_incs[1] + m_incs[2] + 1] += w[7];


				return 0;
			}

		
		return a*fpos1+c*fpos2;//a



	//	if (pos[0]<=1 || pos[0]>=m_dims[0]-2 || pos[1]<=1 || pos[1]>=m_dims[1]-2 || pos[2]<=1 || pos[2]>=m_dims[2]-2) return 0.0f;

	//	float a, b, c, d;
	//	int ipos[3];
	//	float fpos1, fpos2;

	//	float* ptr = m_data;// 这个data值要每一次projection后更新。 SinglesCALARvOLUME data


	//	ipos[0] = (int)pos[0]; ipos[1] = (int)pos[1]; ipos[2] = (int)pos[2];
	//	//
	//	float t0 = pos[0] - ipos[0];
	//	float t1 = pos[1] - ipos[1];
	//	float t2 = pos[2] - ipos[2];
	//	float M_t0 = 1.0f - t0;
	//	float M_t1 = 1.0f - t1;
	//	float M_t2 = 1.0f - t2;
	//	float w[8] = { M_t0*M_t1*M_t2, t0*M_t1*M_t2, M_t0*t1*M_t2, t0*t1*M_t2, M_t0*M_t1*t2, t0*M_t1*t2, M_t0*t1*t2, t0*t1*t2 };

	//	//7-6-5-4-3-2-1-0 weight
	//	/*float _p=(float)t0*t1*t2+M_t0*t1*t2+t0*t2*M_t1+M_t0*M_t1*t2+t0*t1*M_t2+M_t0*t1*M_t2+t0*M_t1*M_t2+M_t0*M_t1*M_t2;
	//	if(!m_IsBackProjectionMode)
	//	SumWeight_Ray+=_p;*/
	//	float aaaa=w[0]+w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7];
	///*	if (aaaa>1.00001||a<0.99999)
	//	{
	//		cout<<"weight:	"<<(double)aaaa<<endl;
	//	}*/
	//	
	//	if (!m_IsBackProjectionMode)
	//		SumWeight_Ray+=aaaa;
	//		//SumWeight_Ray += 1.0f;
	//	//modificated by GM 11/23/2015 speed up
	//	//	if(!m_IsBackProjectionMode)
	//	//	SumWeight_Ray+=1.0f;
	//	//cout<<_p<<endl;
	//	//cout<<SumWeight_Ray<<endl;
	//	int m_loc = ipos[0] + ipos[1] * m_incs[1] + ipos[2] * m_incs[2];
	//	ptr += m_loc;
	//	////int t=ipos[0]+ipos[1]*m_incs[1]+ipos[2]*m_incs[2];
	//	//fpos2=pos[0]-ipos[0];
	//	//fpos1=1.0f-fpos2;

	//	//a=ptr[0]*fpos1+ptr[1]*fpos2;
	//	//b=ptr[m_incs[1]]*fpos1+ptr[m_incs[1]+1]*fpos2;
	//	//c=ptr[m_incs[2]]*fpos1+ptr[m_incs[2]+1]*fpos2;
	//	//d=ptr[m_incs[1]+m_incs[2]]*fpos1+ptr[m_incs[1]+m_incs[2]+1]*fpos2;  

	//	//fpos2=pos[1]-ipos[1];
	//	//fpos1=1.0f-fpos2;

	//	//a=a*fpos1+b*fpos2;//[(ptr[0]*fpos1+ptr[1]*fpos2)  *  (pos[1]-ipos[1])  +  (ptr[m_incs[1]]*(1.0f-fpos2)+ptr[m_incs[1]+1]*fpos2)  *  (pos[1]-ipos[1])]* (1.0f-fpos2)
	//	//c=c*fpos1+d*fpos2;//[(ptr[m_incs[2]]*fpos1+ptr[m_incs[2]+1]*fpos2) * (1.0f-fpos2) + (ptr[m_incs[1]+m_incs[2]]*fpos1+ptr[m_incs[1]+m_incs[2]+1]*fpos2)  *  (pos[1]-ipos[1])]  * (pos[2]-ipos[2])

	//	//fpos2=pos[2]-ipos[2];
	//	//fpos1=1.0f-fpos2;



	//	float gw = ptr[0] * w[0] + ptr[1] * w[1] + ptr[m_incs[1]] * w[2] + ptr[m_incs[1] + 1] * w[3] +
	//		ptr[m_incs[2]] * w[4] + ptr[m_incs[2] + 1] * w[5] + ptr[m_incs[1] + m_incs[2]] * w[6] + ptr[m_incs[1] + m_incs[2] + 1] * w[7];

	//	if (m_IsBackProjectionMode)
	//	{
	//		float* pCtrb = pContributions;
	//		//cout<<m_incs[0]<<"		"<<m_incs[1]<<"		"<<"	"<<m_incs[2]<<endl;
	//		//if(pCorrection[indexOfimage]>2.0)
	//		//cout<<"pCorrection: "<<  pCorrection[indexOfimage]<<endl;
	//		pCtrb += m_loc;
	//		pCtrb[0] += (float)(pCorrection[indexOfimage] * w[0]);
	//		pCtrb[1] += (float)(pCorrection[indexOfimage] * w[1]);
	//		pCtrb[m_incs[1]] += (float)(pCorrection[indexOfimage] * w[2]);
	//		pCtrb[m_incs[1] + 1] += (float)(pCorrection[indexOfimage] * w[3]);
	//		pCtrb[m_incs[2]] += (float)(pCorrection[indexOfimage] * w[4]);
	//		pCtrb[m_incs[2] + 1] += (float)(pCorrection[indexOfimage] * w[5]);
	//		pCtrb[m_incs[1] + m_incs[2]] += (float)(pCorrection[indexOfimage] * w[6]);
	//		pCtrb[m_incs[1] + m_incs[2] + 1] += (float)(pCorrection[indexOfimage] * w[7]);
	//		if (M_t0*M_t1*M_t2>1.0f || t0*M_t1*M_t2>1.0 || M_t0*t1*M_t2>1.0 || t0*t1*M_t2>1.0 || M_t0*M_t1*t2>1.0 || t0*t2*M_t1>1.0 || M_t0*t1*t2>1.0 || t0*t1*t2>1.0)
	//			cout << M_t0*M_t1*M_t2 << " " << t0*M_t1*M_t2 << "	" << M_t0*t1*M_t2 << "" << t0*t1*M_t2 << " " << M_t0*M_t1*t2 << " " << t0*t2*M_t1 << " " << M_t0*t1*t2 << " " << t0*t1*t2 << endl;
	//		
	//		float *pweight = pWeightForVoxels;
	//		pweight += m_loc;

	//		pweight[0] += w[0];
	//		pweight[1] += w[1];
	//		pweight[m_incs[1]] += w[2];
	//		pweight[m_incs[1] + 1] += w[3];
	//		pweight[m_incs[2]] += w[4];
	//		pweight[m_incs[2] + 1] += w[5];
	//		pweight[m_incs[1] + m_incs[2]] += w[6];
	//		pweight[m_incs[1] + m_incs[2] + 1] += w[7];


	//		return 0;
	//	}
	//	//return a*fpos1+c*fpos2;//a

	//	return gw;




	}


	virtual float GetScalarValue_Kaiser(Vector4 pos, float theta)
	{

		//	Vector4 norm_t2(came_test);
		//	norm_t2.Normalize();
		//	//cout<<"Cameradirection:			"<<came_test[0]<<"   "<<came_test[1]<<"   "<<came_test[2]<<"   "<<came_test[3]<<endl;
		//	cout<<came_test*(plane_norm)/CameraDirection.Length()<<"   "<<(float) abs(norm_t2.ele[2]/norm_t2.Length())<<endl;
		// 看一下PRro工程里面的描述。
		return 0;
	}

	virtual Vector4 GetNormal(Vector4 pos)
	{
		return Vector4(0,0,0,1);
	
	}

};

// static volume access
static VolumeAccess *va;


// core function, get the hray value from data
static inline float GetScalarValue(int x,int y,int z,float &SumWeight_Ray,int indexOfimage)
//static inline float GetScalarValue(int x,int y,int z)
{
//	if (va) return va->GetScalarValue(x,y,z); else return 0.0f;
		if (va) return va->GetScalarValue(x,y,z,SumWeight_Ray, indexOfimage); else return 0.0f;
}


// in the case of interpolation
static inline float GetScalarValue(Vector4 pos,float &SumWeight_Ray,int indexOfimage)
//static inline float GetScalarValue(Vector4 pos)
{
	if (va) return va->GetScalarValue(pos,SumWeight_Ray, indexOfimage); else return 0.0f;
}

// in the case of interpolation

 static float GetScalarValue_Kaiser(Vector4 pos, float theta)
{

	if (va) return va->GetScalarValue_Kaiser(pos,theta); else return 0.0f;
}

// get normal
static inline Vector4 GetNormal(Vector4 pos)
{
	if (va) return va->GetNormal(pos); else return Vector4();
}

// Preparation
//static bool SingleScalarVolumePrepare(VolumeData *volume,TypicalRayProducer *rayproducer, int frame, int level)
static bool SingleScalarVolumePrepare(VolumeData *volume, TypicalRayProducer *rayproducer, int frame, bool isWarpped, int ProjId, float rtime, int flevel)
{
	if (!volume) return false;
	//VolumeData::DataType type=volume->GetDataType();
//	float type=volume->GetDataType();
	va=0;
	//if(va) 

	float *data;
	
	int dims[3];
	float spacings[3];
	bool IsBackProMode;
	float *top;
	float *down;
	float *correction;
	//data=volume->GetData();
	//top=volume->GetContributions();
	//down=volume->GetTotalWeight();
	top=rayproducer->GetContributions();
	down=rayproducer->GetTotalWeight();
	volume->GetDimensions(dims);
	volume->GetSpacings(spacings);
	IsBackProMode=volume->GetProjectionMode();
	//correction=volume->GetCorrection();
	correction=rayproducer->GetCorrection();
	image im = volume->m_Frameslist[flevel][frame].fVolume;
	if (isWarpped)
	{
		//int rFrameidx = m_Data->m_ProjIndex[projIter].rightframe;
		float rTime =  rtime;
		float aabb[] = { 0,dims[0],0,dims[1]  ,0,dims[2] };
		int nfames = volume->m_Frameslist[flevel].size();
if(rTime<0.0f)
		im = advection::advect(dims, aabb, volume->m_Frameslist[flevel][max(0, frame - 1)].fFlow, volume->m_Frameslist[flevel][frame].fVolume, rTime);
else
		im = advection::advect(dims, aabb, volume->m_Frameslist[flevel][min(nfames-1, frame)].fFlow, volume->m_Frameslist[flevel][frame].fVolume, rTime);

	}
	//volume->m_Frameslist[frame].
	 //va=new T_VolumeAccess(data,dims,spacings,IsBackProMode,top,down,correction);
	// va = new T_VolumeAccess(volume->m_Volumelist[frame], dims, spacings, IsBackProMode, top, down, correction);
	 //va = new T_VolumeAccess(volume->m_VolumeLevel[level][frame], dims, spacings, IsBackProMode, top, down, correction);
	// va = new T_VolumeAccess(volume->m_Frameslist[frame].fVolume, dims, spacings, IsBackProMode, top, down, correction);
	 va = new T_VolumeAccess(im, dims, spacings, IsBackProMode, top, down, correction);
	 


	/* idata = volume->m_Volumelist[frame];
	 pContributions = (float*)top;
	 pWeightForVoxels = (float*)down;
	 pCorrection = (float*)correction;
	 m_IsBackProjectionMode = IsBackProMode;
	 m_dims[0] = dims[0];m_dims[1] = dims[1];m_dims[2] = dims[2];
	 m_spacings[0] = spacings[0];m_spacings[1] = spacings[1];m_spacings[2] = spacings[2];

	 m_incs[0] = 1;
	 m_incs[1] = m_dims[0];
	 m_incs[2] = m_dims[0] * m_dims[1];*/


	//if (type.isFloat) va=new T_VolumeAccess<float>(data,dims,spacings);
	//else if (type.isSigned)
	//{
	//	if (type.bitsPerSample==8) va=new T_VolumeAccess<char>(data,dims,spacings);
	//	else if (type.bitsPerSample==16) va=new T_VolumeAccess<short>(data,dims,spacings);
	//	else if (type.bitsPerSample==32) va=new T_VolumeAccess<int>(data,dims,spacings);
	//}
	//else
	//{
	//	if (type.bitsPerSample==8) va=new T_VolumeAccess<unsigned char>(data,dims,spacings);
	//	else if (type.bitsPerSample==16) va=new T_VolumeAccess<unsigned short>(data,dims,spacings);
	//	else if (type.bitsPerSample==32) va=new T_VolumeAccess<unsigned int>(data,dims,spacings);
	//}

	return true;
}

// recycling
static void SingleScalarVolumeClear()
{
	if (va) delete va;
}








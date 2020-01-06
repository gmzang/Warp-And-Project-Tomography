// VolumeData.cpp

#include "VolumeData.h"
#include <malloc.h>
#include <fstream>
#include <iostream>
using namespace  std;

typedef std::vector<Frame_item> framelist;


VolumeData::VolumeData()
{
	
	m_Dims[0]=128;
	m_Dims[1]=128;
	m_Dims[2]=128;
	m_Spacings[0]=m_Spacings[1]=m_Spacings[2]=1.0f;


	int _index=0 ;
	
	m_IsBackProjectMode=false;
	
}


 VolumeData::VolumeData(int d1,int d2, int d3, float sp1,float sp2,float sp3,float *data)
{
	int d[]={d1,d2,d3};

	this->SetDimensions(d);
	float sp[]={sp1,sp2,sp3};
	this->SetSpacings(sp);
	this->_allocate();
//	this->SetSpacings[]

}


 //void VolumeData::setParas(int vx,int vy, int vz,double vs)

 void VolumeData::setParas(int *vols, double vs,int nframes, int totalproj)
 {
	 m_volX= vols[0];
	 m_volY= vols[1];
	 m_volZ= vols[2];
	 Vol_Spacing=vs;
	 m_nframes = nframes;
	 m_totalProjs = totalproj;

 }
 
VolumeData::~VolumeData()
{
	Clean();
}

void VolumeData::Clean()
{
	/*if (m_Data)
	{
		free(m_Data);
		m_Data=0;
	}*/

	
}


void VolumeData::_allocate()
{

}



void VolumeData::GetDimensions(int *dims) const
{
	dims[0]=m_Dims[0]; dims[1]=m_Dims[1]; dims[2]=m_Dims[2];
}

void VolumeData::SetDimensions( int *dims)
{
	m_Dims[0]=dims[0]; m_Dims[1]=dims[1]; m_Dims[2]=dims[2];
}

void VolumeData::GetSpacings(float *spacings) const
{
	spacings[0]=m_Spacings[0]; spacings[1]=m_Spacings[1]; spacings[2]=m_Spacings[2];
}

void VolumeData::SetSpacings(const float *spacings)
{
	m_Spacings[0]=spacings[0]; m_Spacings[1]=spacings[1]; m_Spacings[2]=spacings[2];
}




CImg<float> VolumeData::GetData(int iframe)
{
	return m_Frameslist[0][0].fVolume;

}


void  VolumeData::SetData()
{
	//m_Data = _data;
	//return true;
	//m_Volumelist = _vol_list;

}


//void VolumeData::ReadStructuredFile(const char* filename)
void VolumeData::ReadStructuredFile()
{

	/*FILE *fp;
	fp=fopen(filename,"rb");*/


	std::vector<Frame_item> m_List[5];
	m_Dims[0]=m_volX;
	m_Dims[1]=m_volY;
	m_Dims[2]=m_volZ;

	m_Spacings[0]=m_Spacings[1]=m_Spacings[2]=Vol_Spacing;  //0.3751f


	float baseNo = 3.0f;
	int nLevel = 2;
	int maxFrames = 130;
	size_item  szitem;

	szitem.volsize[0] = m_Dims[0] ;
	szitem.volsize[1] = m_Dims[1] ;
	szitem.volsize[2] = m_Dims[2] ;
	cout << "szitem.volsize: " << szitem.volsize[0] << " " << szitem.volsize[1] << " " << szitem.volsize[2] << endl;
	m_Volumesize.push_back(szitem);
	
	//float timeunit = 1.0f;
	//float timeunit = (float)((float)m_totalProjs / (float)m_nframes);
	

//	float fistslot = (timeunit - 1) / 2.0f;
	//cout << "timeunit: " << m_proNo<<" "<< timeunit<< endl;
	//cout << "Time unit:	" << timeunit << " " << endl;
	//cout << "Time unit:	" << fistslot << " " << endl;
	char basename[1024];
	sprintf(basename, "input");

	int factor = 1;
	for (int  ilevel = 1; ilevel <= nLevel; ilevel++)
	{

		int nf = min(m_nframes*factor, maxFrames);
		m_fInLevel[ilevel - 1] = nf;
		float timeunit = (float)((float)m_totalProjs / (float)nf);
		m_tUnitInLevel[ilevel - 1] = timeunit;
		//baseNo =timeunit;

		baseNo = max(min(timeunit,6.0f),1.0f)*1.0/float(factor);
		baseNo = max(baseNo, 1.0f);
		baseNo = 3;
		m_base = baseNo;
		cout << nf << endl;
		//for (int i = 0; i < m_nframes; i++)
		for (int  i = 0; i < nf; i++)
		{
			
			Frame_item item;
			//int sz[] = { m_Dims[0] * downsample_eta ,m_Dims[1] * downsample_eta,m_Dims[2] * downsample_eta };
			int sz[] = { m_Dims[0]  ,m_Dims[1] ,m_Dims[2]  };
			item.SetVolSize(sz);
			item.fVolume = image(sz[0], sz[1], sz[2], 1, 0.0f);
			//item.fVolume = image(temp);

			//item.fFlow = image(sz[0], sz[1], sz[2], 3, 0.0f);
			
			// some issues for relation between fTime and item.fStartproj
			item.fTime = timeunit * i;

			//item.fStartproj = max(0, item.fTime - (int)((float)item.fTime*0.5 / (float)baseNo));
			////item.fEndproj = min(item.fStartproj + (int)baseNo, m_proNo);
			//item.fEndproj = min(item.fTime + (int)((float)item.fTime*0.5 / (float)baseNo), m_totalProjs);
	

			item.fStartproj = max(0, (int)(item.fTime - baseNo  ));
			//item.fStartproj = max(0, (int)(item.fTime - (baseNo/2.0))+1);
			//item.fEndproj = min(item.fStartproj + (int)baseNo, m_proNo);
			//item.fEndproj = min((int)(item.fStartproj+ baseNo-1), m_totalProjs - 1);
			item.fEndproj = min((int)(item.fTime + baseNo ), m_totalProjs - 1);
			//item.fEndproj = min((int)(item.fTime + (baseNo / 2.0)), m_totalProjs-1);
			item.fprojRange = item.fEndproj - item.fStartproj;
	
			cout <<"item: "<<i +1<< "  item.fStartproj:	" << item.fStartproj << " <--- " << item.fTime << " ---> " << item.fEndproj << ", Range:" << item.fprojRange << endl;
			m_Frameslist[ilevel-1].push_back(item);
			//item.fsize = sz;
		}
		m_Frameslist[ilevel - 1].push_back(m_Frameslist[ilevel - 1][nf-1]);
		m_Frameslist[ilevel - 1].push_back(m_Frameslist[ilevel - 1][nf - 1]);
		//system("pause");
		//m_List[2][2] = m_List[1][2];

		factor *= 2;

	}






	//print basic information
	std::cout <<"ReadStructuredFile:	"<<m_Dims[0] <<"  "<<m_Dims[1]<<"   "<<m_Dims[2]<< " frames:  " << m_nframes<<endl;
	std::cout << m_Spacings[0] <<"  "<<m_Spacings[1]<<"   "<<m_Spacings[2]<<endl;
	/*for (int j = 0; j < m_Dims[0] * m_Dims[1] * m_Dims[2]; j++)
	{
		if (m_Data[j] >1.1) cout << m_Data[j] << endl; 
	}*/
	//std::cout << " voxel values:" <<m_Data[1] << " " << m_Data[1000] << " " << m_Data[10000] << " " << m_Data[100000] << endl;
	//std::cout << m_Dims[0]*m_Dims[1]*m_Dims[2]<<endl;
	std::cout<<"==========Finish the input volume information========="<<endl;

}


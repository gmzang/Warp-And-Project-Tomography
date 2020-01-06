//#include "stdafx.h"
#include "ImageIO.h"
#include <iostream>
using namespace  std;
/*------------------------------------------------
*函数名称：
*	readBmp()
*
*函数参数：
*	char *bmpName -文件名字及路径
*	int *bmpWidth -图像宽度指针变量
*	int *bmpHeight -图像高度指针变量
*
*返回值：
*	unsigned char * -有效图像数据指针
*
*说明：
*	（1）给定一个图像文件名及其路径，读图像的位图数据、宽、高、颜色表及每像素
*位数等数据进内存。
*	（2）该读程序仅针对灰度图像（biBitCount=8）格式
------------------------------------------------*/
BYTE* Read8bitbmp(char *bmpName,int *bmpWidth,int *bmpHeight)
{
	//二进制读方式打开指定的图像文件
	FILE *fp=fopen(bmpName,"rb");
	if(fp==0) return 0;

	//跳过位图文件头结构BITMAPFILEHEADER
	fseek(fp, sizeof(BITMAPFILEHEADER),0);	//moves the file pointer to a specified location

	//定义位图信息头结构变量，读取位图信息头进内存，存放在变量head中
	BITMAPINFOHEADER head;
	fread(&head, sizeof(BITMAPINFOHEADER), 1,fp); 

	//获取图像宽、高等信息
	*bmpWidth = head.biWidth;	//位图的宽度，以像素为单位
	*bmpHeight = head.biHeight;	//位图的高度，以像素为单位

	//定义变量，计算图像每行像素所占的字节数（必须是4的倍数）
	int biBitCount=8;	//8位灰度图像
	int lineByte=((*bmpWidth)*biBitCount/8+3)/4*4;

	//灰度图像有颜色表，且颜色表表项为256	
	RGBQUAD *pColorTable=new RGBQUAD[256];	//申请颜色表所需要的空间
	fread(pColorTable,sizeof(RGBQUAD),256,fp);	//读颜色表进内存

	//申请位图数据所需要的空间，读位图数据进内存
	unsigned char *pBmpBuf=new unsigned char[lineByte*(*bmpHeight)];
	fread(pBmpBuf,1,lineByte*(*bmpHeight),fp);
	
	//关闭文件
	fclose(fp);

	//提取有效数据，并返回
	if (*bmpWidth<lineByte)
	{
		unsigned char *pBmpBufData=new unsigned char[(*bmpWidth)*(*bmpHeight)];
		for (int i=0;i<*bmpHeight;i++)
		{
			for (int j=0;j<*bmpWidth;j++)
			{
				pBmpBufData[i*(*bmpWidth)+j]=*(pBmpBuf+i*lineByte+j);	//第一种赋值方法
				//*(pBmpBufData+i*(*bmpWidth)+j)=*(pBmpBuf+i*lineByte+j);	//第二种赋值方法
			}
		}
		return pBmpBufData;	//返回提取的有效数据
	}
	else
	{
		return pBmpBuf;	//返回有效数据
	}
}

/*------------------------------------------------
*函数名称：
*	saveBmp()
*
*函数参数：
*	unsigned char *imgBuf-待存盘的位图数据
*   <span style="WHITE-SPACE: pre">	</span>int width-以像素为单位待存盘位图的宽
*   <span style="WHITE-SPACE: pre">	</span>int height-以像素为单位待存盘位图高
*	char *bmpName-文件名字及路径
*
*返回值：
*   0为失败，1为成功
*
*说明：
*	1.给定一个图像的位图数据、宽、高等信息，将其写到指定文件中。
*	2.该读程序仅针对灰度图像（biBitCount=8）格式
------------------------------------------------*/
bool Writebitmapbmp(unsigned char *imgBuf, int width, int height, char *bmpName) 
{
	//如果位图数据指针为0，则没有数据传入，函数返回
	if(!imgBuf) return 0;

	int biBitCount=8;	//每个像素所占的位数(bit)

	//颜色表大小，以字节为单位，灰度图像颜色表为1024字节
	int colorTablesize=1024;

	//待存储图像数据每行字节数为4的倍数
	int lineByte=(width * biBitCount/8+3)/4*4;

	//以二进制写的方式打开文件
	FILE *fp=fopen(bmpName,"wb");
	if(fp==0) return 0;

	//申请位图文件头结构变量，填写文件头信息
	BITMAPFILEHEADER fileHead;
	fileHead.bfType=0x4D42;	//bmp类型

	//bfSize是图像文件4个组成部分之和
	fileHead.bfSize=sizeof(BITMAPFILEHEADER)+sizeof(BITMAPINFOHEADER)+colorTablesize+lineByte*height;
	fileHead.bfReserved1=0;
	fileHead.bfReserved2=0;

	//bfOffBits是图像文件前3个部分所需空间之和
	fileHead.bfOffBits=54+colorTablesize;

	//写文件头进文件
	fwrite(&fileHead, sizeof(BITMAPFILEHEADER),1, fp);

	//申请位图信息头结构变量，填写信息头信息
	BITMAPINFOHEADER head;
	head.biSize=40;	//本结构的长度
	head.biWidth=width;	//位图的宽度，以像素为单位
	head.biHeight=height;	//位图的宽度，以像素为单位
	head.biPlanes=1;	//目标设备的级别，必须是1
	head.biBitCount=8;	//每个像素所占的位数(bit)，其值必须为1（黑白图像），4（16色图），8（256色），24（真彩色图）
	head.biCompression=BI_RGB;	//位图压缩类型，有效的值为BI_RGB（未经压缩）、BI_RLE4，BI_RLE8，BI_BITFIEDS（均为Windows定义常量）。
	head.biSizeImage=lineByte*height;	//实际的位图数据占用的字节数
	head.biXPelsPerMeter=0;	//指定目标数据的水平分辨率，单位是像素/米。
	head.biYPelsPerMeter=0;	////指定目标数据的垂直分辨率，单位是像素/米。
	head.biClrUsed=0;	//位图实际用到的颜色数，如果该值为零，则用到的颜色数为2的biBitCount次幂
	head.biClrImportant=0;	//位图显示过程中重要的颜色数，如果该值为零，则认为所有的颜色都是重要的。

	//写位图信息头进内存
	fwrite(&head, sizeof(BITMAPINFOHEADER),1, fp);

	//申请颜色表所需要的空间，写颜色表进文件
	RGBQUAD *pColorTable=new RGBQUAD[256];
	for (int i=0;i<256;i++)
	{
		pColorTable[i].rgbRed=i;
		pColorTable[i].rgbGreen=i;
		pColorTable[i].rgbBlue=i;
		pColorTable[i].rgbReserved=0;
	}
	fwrite(pColorTable,sizeof(RGBQUAD),256,fp);

	//判断位图数据宽度是否正确，并写数据进入BMP
	if (width<lineByte)	//如果有效数据宽度小于BMP格式要求宽度
	{
		unsigned char *imgBufBMP=new unsigned char [lineByte*height];
		//将有效数据赋给用于写BMP的数据空间
		for (int i=0;i<height;i++)
		{
			for (int j=0;j<width;j++)
			{
				*(imgBufBMP+i*lineByte+j)=*(imgBuf+i*width+j);	
			}
		}
		//将大于有效数据宽度的部分补0
		for (int i=0;i<height;i++)
		{
			for (int j=width;j<lineByte;j++)
			{
				*(imgBufBMP+i*lineByte+j)=0;	
			}
		}
		//写BMP格式要求的图数据写进文件
		fwrite(imgBufBMP, height*lineByte, 1, fp);
	}
	else
	{
		//写有效图数据写进文件
		fwrite(imgBuf, height*lineByte, 1, fp);
	}

	//关闭文件
	fclose(fp); 

	return 1;
}



////#include "stdafx.h"
////#include "stdafx.h" 
//#include "ImageIO.h"
// 
//#include "Windows.h"
//bool ImageIO::readBmp(char *bmpName)
//{
//	//二进制读方式打开指定的图像文件
//	FILE *fp=fopen(bmpName,"rb");
//	if(fp==0) return 0;
//
//
//	//跳过位图文件头结构BITMAPFILEHEADER
//	fseek(fp, sizeof(BITMAPFILEHEADER),0);
//
//	//定义位图信息头结构变量，读取位图信息头进内存，
//	//存放在变量head中
//	BITMAPINFOHEADER head;  
//	fread(&head, sizeof(BITMAPINFOHEADER), 1,fp); 
//	//获取图像宽、高、每像素所占位数等信息
//	bmpWidth = head.biWidth;
//	bmpHeight = head.biHeight;
//	biBitCount = head.biBitCount;
//	//定义变量，计算图像每行像素所占的字节数（必须是4的倍数）
//	int lineByte=(bmpWidth * biBitCount/8+3)/4*4;
//	//灰度图像有颜色表，且颜色表表项为256
//	if(biBitCount==8){
//		//申请颜色表所需要的空间，读颜色表进内存
//		pColorTable=new RGBQUAD[256];
//		fread(pColorTable,sizeof(RGBQUAD),256,fp);
//	}
//	//申请位图数据所需要的空间，读位图数据进内存
//	pBmpBuf=new unsigned char[lineByte * bmpHeight];
//	fread(pBmpBuf,1,lineByte * bmpHeight,fp);
//	//关闭文件
//	fclose(fp);
//	return 1;
//}
//
//bool ImageIO::saveBmp(char *bmpName, unsigned char
//	*imgBuf, int width, int height, 
//	int biBitCount, RGBQUAD *pColorTable)
//{
//	//如果位图数据指针为0，则没有数据传入，函数返回
//	if(!imgBuf)
//		return 0;
//	//颜色表大小，以字节为单位，灰度图像颜色表
//	//	为1024字节，彩色图像颜色表大小为0
//	int colorTablesize=0;
//	if(biBitCount==8)
//		colorTablesize=1024;
//	//待存储图像数据每行字节数为4的倍数
//	int lineByte=(width * biBitCount/8+3)/4*4;
//	//以二进制写的方式打开文件
//	FILE *fp=fopen(bmpName,"wb");
//	if(fp==0) return 0;
//	//申请位图文件头结构变量，填写文件头信息
//	BITMAPFILEHEADER fileHead;
//	fileHead.bfType = 0x4D42;//bmp类型
//	//bfSize是图像文件4个组成部分之和
//	fileHead.bfSize= sizeof(BITMAPFILEHEADER)
//		+ sizeof(BITMAPINFOHEADER)
//		+ colorTablesize + lineByte*height;
//	fileHead.bfReserved1 = 0;
//	fileHead.bfReserved2 = 0;
//	//bfOffBits是图像文件前3个部分所需空间之和
//	fileHead.bfOffBits=54+colorTablesize;
//	//写文件头进文件
//	fwrite(&fileHead, sizeof(BITMAPFILEHEADER),1, fp);
//	//申请位图信息头结构变量，填写信息头信息
//	BITMAPINFOHEADER head; 
//	head.biBitCount=biBitCount;
//	head.biClrImportant=0;
//	head.biClrUsed=0;
//	head.biCompression=0;
//	head.biHeight=height;
//	head.biPlanes=1;
//	head.biSize=40;
//	head.biSizeImage=lineByte*height;
//	head.biWidth=width;
//	head.biXPelsPerMeter=0;
//	head.biYPelsPerMeter=0;
//	//写位图信息头进内存
//	fwrite(&head, sizeof(BITMAPINFOHEADER),1, fp);
//	//如果灰度图像，有颜色表，写入文件 
//	if(biBitCount==8)
//		fwrite(pColorTable, sizeof(RGBQUAD),256, fp);
//	//写位图数据进文件
//	fwrite(imgBuf, height*lineByte, 1, fp);
//	//关闭文件
//	fclose(fp);
//	return 1;
//}
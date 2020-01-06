#ifndef _ReadWriteBmp_H
#define _ReadWriteBmp_H
#include <windows.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

BYTE *Read8bitbmp(char *bmpName,int *bmpWidth,int *bmpHeight);	//read bmp file 
bool Writebitmapbmp(unsigned char *imgBuf, int width, int height, char *bmpName);	// write bmp file

#endif



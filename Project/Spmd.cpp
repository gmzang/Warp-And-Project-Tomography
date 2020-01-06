//// SPMD.cpp
//#include "SPMD.h"
//#include "pthread.h"
//#include "NumberOfProcessors.h"
//
//struct ThreadInfo
//{
//	SPMDThread entrance;
//	void *userData;
//	int NumberOfThreads;
//	int ThreadID;
//	pthread_t SystemThreadID;
//};
//
//static void* _thread(void* info)
//{
//	ThreadInfo* tInfo=(ThreadInfo*)info;
//	tInfo->entrance(tInfo->NumberOfThreads,tInfo->ThreadID,tInfo->userData);
//	return 0;
//}
//
//void LaunchSPMD(SPMDThread entrance,void* userData)
//{
//	int nop=GetNumberOfProcessors();
//	int i;
//
//	ThreadInfo* tInfo=new ThreadInfo[nop];
//
//	for (i=1;i<nop;i++)
//	{
//		tInfo[i].entrance=entrance;
//		tInfo[i].userData=userData;
//		tInfo[i].NumberOfThreads=nop;
//		tInfo[i].ThreadID=i;
//		int threadError = pthread_create(&tInfo[i].SystemThreadID,NULL,_thread,&tInfo[i]);
//	}
//
//	entrance(nop,0,userData);
//
//	for (i=1;i<nop;i++)
//		pthread_join(tInfo[i].SystemThreadID,NULL);
//
//	delete[] tInfo;
//}
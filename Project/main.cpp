#include "GLUTSceneContainer.h" // 

#include "WaP_Tomo.h"
#include "VolumeData.h" // 数据
#include <argtable2.h>
//#include<iostream>
//using namespace std;
extern "C" { FILE __iob_func[3] = { *stdin,*stdout,*stderr }; }

int main(int argc,char** argv)
{



	struct arg_dbl *d1=arg_dbl0("s","Sigma",        "<double>"     ,     "Mu in PSART framework(default as 0.2)");
	struct arg_dbl *d2=arg_dbl0("t","Tau",  "<double>"       ,         "Lambda in PSART framework(default as 0.2)");
	//  struct arg_lit  *list    = arg_lit0("lL",NULL,                      "list files");
	//struct arg_lit  *recurse = arg_lit0("R",NULL,                       "recurse through subdirectories");
	struct arg_dbl *d3=arg_dbl0("l","Lambda",  "<double>"       ,         "trade-off weight for data and prior ");
	struct arg_dbl *alp=arg_dbl0("a","Alpha",  "<double>"       ,         "weight for sart algorithm to control the convergence rate default 0.1 ");
	
	struct arg_file *outfile = arg_file0("o",NULL,"<output>",           "output file (default is \"hello.txt\")");
	//struct arg_int *proNo=arg_int0("n","imageNo",  "<int>"       ,         "Number of projection images default 360");

	//struct arg_int *proNo=arg_int0("n","imageNo",  "<int>"       ,         "Number of Cone beam projection images");
	struct arg_int *nframe = arg_int0("v", "nframes", "<int>", "Number of frames");
	struct arg_int *nrounds = arg_int0("r", "nframes", "<int>", "Number of rounds for each frames");



	struct arg_int  *vols = arg_intn("f", "XYZ", NULL, 1, 3, "takes an integer value (defaults to 9)");

	//struct arg_int  *voltest= arg_int0("e", "fot", NULL,  "takes an integer value (defaults to 9)");

	struct arg_int *imgWH = arg_intn("u", "imgUV", NULL,1,2, "w and h of projection image default:512 ");

	//struct arg_int *imgW=arg_int0("u","imageW",  "<int>"       ,         "width of projection image default:512 ");
	//struct arg_int *imgH=arg_int0("v","imageH",  "<int>"       ,         "height of projection image  default 512");

	//struct arg_int *volX=arg_int0("x","vdX",  "<int>"       ,         "volume size (default as 256)");
	//struct arg_int *volY=arg_int0("y","vdY",  "<int>"       ,         "volume size (default as 256)");
	//struct arg_int *volZ=arg_int0("z","vdZ",  "<int>"       ,         "volume size (default as 256)");
	struct arg_dbl *vs=arg_dbl0("c","voxelspacing",  "<double>"       ,         "voxel spacing (default as 1)");
	struct arg_dbl *ds=arg_dbl0("d","detector spacing",  "<double>"       ,         "detector spacing (default as 1)");
	//struct arg_file *projfile = arg_file0("i",NULL,"<proj>",           "Input Dataset file(.mha\.tif\.bmp and other raw file are supported)");
	struct arg_file *projfiles = arg_filen("i","projsvolume", NULL,1,20, "Input Dataset file(.mha\.tif\.bmp\.hdr and other raw file are supported)");

	struct arg_dbl *sdd=arg_dbl0("j","sdd",  "<double>"       ,         " Source to Detector Distance (default as 1000mm)");
	//struct arg_dbl *ox=arg_dbl0("k","xx",  "<double>"       ,         " X offset on projection image  (default as 0)");
	//struct arg_dbl *oy=arg_dbl0("m","yy",  "<double>"       ,         " y offset on projection image (default as 0)");
	struct arg_dbl *sid=arg_dbl0("g","sid",  "<double>"       ,         " Source to Iso-Object Distance  (default as 600mm)");
	//struct arg_dbl *cone_beam_angle=arg_dbl0("r","arc",  "<double>"       ,         " Cone beam angle ,default as 25 degree");

	struct arg_dbl *offsetXYZ = arg_dbln(NULL, "oxyz", NULL, 1, 3, " offset x y z of the center of volume (default as 0)");
	struct arg_dbl *startdegree = arg_dbln(NULL, "sdg", NULL, 1, 20, " start degree for each proj sequence");
//	struct arg_dbl *RoundsStartdegree = arg_dbln(NULL, "rsdg", NULL, 1, 150, " start degree for each proj sequence");
	/*struct arg_dbl *offsetx=arg_dbl0(NULL,"ox",  "<double>"       ,         " offset x of the center of volume (default as 0)");
	struct arg_dbl *offsety=arg_dbl0(NULL,"oy",  "<double>"       ,         " offset y of the center of volume (default as 0)");
	struct arg_dbl *offsetz=arg_dbl0(NULL,"oz",  "<double>"       ,         " offset z of the center of volume (default as 0)");
*/

	struct arg_int *AlgoIter=arg_int0("b","AlgorithmIter",  "<int>"       ,         "Number of algorithm iterations, default as 20");

	struct arg_int *SartIter=arg_int0("p","SartIter",  "<int>"       ,         "Number of SART nested iterations, default as 1");

	// add for prior options
	struct arg_str  *prior    = arg_str0(NULL,"prior", "{STV,SAD,ATV}", "specify the prior you are applying from {STV,SAD,ATV} , default as STV");
   struct arg_str  *bp    = arg_str0(NULL,"bp", "{Voxelbased,Raybased}", "specify the backprojection methods you are applying from {Voxelbased,Raybased}, default as voxelbased");
   


	struct arg_lit  *help    = arg_lit0("h","help",                    "print this help and exit");
	struct arg_lit  *version = arg_lit0(NULL,"version",                 "print version information and exit");

	struct arg_end  *end     = arg_end(20);


	void* argtable[] = {d1,d2,d3,alp,outfile,projfiles,vols,imgWH,vs,ds,sid,sdd,offsetXYZ,startdegree,nframe,nrounds,AlgoIter,SartIter,prior, bp,help,version,end};
	//void* argtable[] = { d1,d2,d3,alp,outfile,projfile,vols,imgW,imgH,vs,ds,sid,sdd,ox,oy,offsetx,offsety,offsetz,proNo,AlgoIter,SartIter,prior, bp,help,version,end };

	const char* progname = "Test";
	int nerrors;
	int exitcode=0;

	/* verify the argtable[] entries were allocated sucessfully */
	if (arg_nullcheck(argtable) != 0)
	{
		/* NULL entries were detected, some allocations m_must have failed */
		printf("%s: insufficient memory\n",progname);
		exitcode=1;
		goto exit;
	}

	/* set any command line default values prior to parsing */
	//repeat->ival[0]=3;
	//projNo->ival[0]=45;
	outfile->filename[0]="Test/test";
	d1->dval[0]=0.3;
	d2->dval[0]=0.3;
	d3->dval[0]=0.1;
	//volX->ival[0]=256;//=volsize->dval[1]=volsize->dval[2]=128;
	//volY->ival[0]=volZ->ival[0]=256;

	//volX->ival[1] = 256;//=volsize->dval[1]=volsize->dval[2]=128;

	//imgH->ival[0]=512;
	//imgW->ival[0]=512;
//	proNo->ival[0]=360;
	AlgoIter->ival[0]=2;
	SartIter->ival[0]=2;
	alp->dval[0]=0.1;
	vs->dval[0]=1.0;
	/*offsetx->dval[0]=0.0;
	offsety->dval[0]=0.0;
	offsetz->dval[0]=0.0;*/
	sdd->dval[0]=200.0f;
	//ox->dval[0]=0.0f;
	//oy->dval[0]=0.0f;
	vols->ival[0] = 63;
	vols->ival[1] = 63;
	vols->ival[2] = 63;
	offsetXYZ->count = 3;
	offsetXYZ->dval[0] = 0.0f;
	offsetXYZ->dval[1] = 0.0f;
	offsetXYZ->dval[2] = 0.0f;
	/*voltest->ival[0] =256 ;
	voltest->ival[1] =256;
	voltest->ival[2] =256;*/
	/*voltest->count = 3;*/
	// add prior
	 prior->sval[0]  = "ATV";     /* --prior={STV,SAD,ATV} */
	 bp->sval[0]  = "Voxelbased";
	//cone_beam_angle->dval[0]=25.0;
	sid->dval[0]=400.0;
	projfiles->filename[0]="proj";
	ds->dval[0]=1.0;
	nframe->ival[0] = 1;
	nrounds->ival[0] = 1;

	//lutR->dval[0]=2.0f;
	//lutfile->filename[0]="L_368_125.bin";
	/* Parse the command line as defined by argtable[] */
	nerrors = arg_parse(argc,argv,argtable);

	/* special case: '--help' takes precedence over error reporting */
	if (help->count > 0)
	{
		printf("Usage: %s", progname);
		arg_print_syntax(stdout,argtable,"\n");
		printf("This program demonstrates the use of the argtable2 library\n");
		printf("for parsing command line arguments. Argtable accepts integers\n");
		printf("in decimal (123), hexadecimal (0xff), octal (0o123) and binary\n");
		printf("(0b101101) formats. Suffixes KB, MB and GB are also accepted.\n");
		arg_print_glossary(stdout,argtable,"  %-25s %s\n");
		exitcode=0;
		goto exit;
	}

	/* special case: '--version' takes precedence error reporting */
	if (version->count > 0)
	{
		printf("'%s' example program for the \"argtable\" command line argument parser.\n",progname);
		printf("September 2003, Stewart Heitmann\n");
		exitcode=0;
		goto exit;
	}

	/* If the parser returned any errors then display them and exit */
	if (nerrors > 0)
	{
		/* Display the error details contained in the arg_end struct.*/
		arg_print_errors(stdout,end,progname);
		printf("Try '%s --help' for more information.\n",progname);
		exitcode=1;
		goto exit;
	}

	/* special case: uname with no command line options induces brief help */
	if (argc==1)
	{
		printf("Try '%s --help' for more information.\n",progname);
		exitcode=0;
		goto exit;
	}



	glutInit(&argc,argv);
	//vols->ival[0] = 181;
	//vols->ival[1] = 250;
	//vols->ival[2] = 181;
	VolumeData* data=new VolumeData;
	//data->setParas(volX->ival[0],volY->ival[0],volZ->ival[0],vs->dval[0]);
	
	// change place 1
	//int totalproj = nrounds->ival[0]* nframe->ival[0]*20; //60
	//int totalproj = nrounds->ival[0] * nframe->ival[0] * 1; //60
	int totalproj = 192; //60

	data->setParas(vols->ival, vs->dval[0], nframe->ival[0], totalproj);

	
	data->ReadStructuredFile(); 
	
	
	
	PSART *scene = new PSART(d1->dval[0], d2->dval[0], d3->dval[0], outfile->filename[0], 
		imgWH->ival, totalproj, AlgoIter->ival[0], SartIter->ival[0], alp->dval[0],
		sid->dval[0], projfiles->filename, offsetXYZ->dval, sdd->dval[0], ds->dval[0], prior->sval[0], 
		bp->sval[0], vols->ival, nframe->ival[0], startdegree->dval,nrounds->ival[0]);

	scene->SetData(data);
	GLUTSceneContainer *mainWindow=GLUTSceneContainer::GetMainWindow(imgWH->ival,ds->dval[0]);

	//GLUTSceneContainer *mainWindow;// = GLUTSceneContainer::GetMainWindow(imgWH->ival, ds->dval[0]);
	mainWindow->SetTitle("WaP Tomo Framework");
	mainWindow->SetScene(scene);

	mainWindow->Create(10, 850);



exit:
	/* deallocate each non-null entry in argtable[] */
	arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
	//system("pause");
	return 0;
}


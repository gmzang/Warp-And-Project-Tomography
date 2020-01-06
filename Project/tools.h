

float dxp( float *data,  int *sizeMat, int i, int j);
float dyp( float *data,  int *sizeMat, int i, int j);
float dxm( float *data,  int *sizeMat, int i, int j);
float dym( float *data,  int *sizeMat, int i, int j);
int index2DtoLinear( int *sizeMat, int i, int j);
int index3DtoLinear( int *sizeMat, int i, int j, int k);
float myAbs(float x);
float myMin(float a, float b);
float myMax(float a, float b);
float myPow2(float x);
float dxc( float *data, float *u,  int *sizeMat, int i, int j);
float dyc( float *data, float *u,  int *sizeMat, int i, int j);
float dxcT( float *data, float *u,  int *sizeMat, int i, int j);
float dycT( float *data, float *u,  int *sizeMat, int i, int j);
float dxp3( float *data,  int *sizeMat, int i, int j, int k);
float dyp3( float *data,  int *sizeMat, int i, int j, int k);
float dzp3( float *data,  int *sizeMat, int i, int j, int k);
float dtp3( float *data,  int *sizeMat, int i, int j, int k);
float dxm3( float *data,  int *sizeMat, int i, int j, int k);
float dym3( float *data,  int *sizeMat, int i, int j, int k);
float dzm3( float *data,  int *sizeMat, int i, int j, int k);
float dtm3( float *data,  int *sizeMat, int i, int j, int k);

float dxc3( float *data,  float *u,  int *sizeMat, int i, int j, int k);
float dyc3( float *data,  float *u,  int *sizeMat, int i, int j, int k);
float dxcT3( float *data,  float *u,  int *sizeMat, int i, int j, int k);
float dycT3( float *data,  float *u,  int *sizeMat, int i, int j, int k);

float dxUpwind( float *data, float *v,  int *sizeMat, int i, int j, int k);
float dyUpwind( float *data, float *v,  int *sizeMat, int i, int j, int k);
float dxUpwindT( float *data, float *v,  int *sizeMat, int i, int j, int k);
float dyUpwindT( float *data, float *v,  int *sizeMat, int i, int j, int k);

int linearTo3Di( int *sizeMat, int index);
int linearTo3Dj( int *sizeMat, int index);
int linearTo3Dk( int *sizeMat, int index);
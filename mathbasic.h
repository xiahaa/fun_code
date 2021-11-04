#ifndef _MATH_BASIC_H
#define _MATH_BASIC_H


struct _CamCoeff
{
	float cx;
	float cy;
	float ux;
	float uy;

	//CamCoeff(){cx = 0.0,cy = 0.0,ux = 0.0,uy = 0.0;};
};

int q2r(const float *pq,   //1x4 quaternions 
		float *pr	   //3x3 rotation matrix	
		);
int r2e(const float *pr,   //3x3 rotation matrix	 
		float *pe	   //3x1 euler angles, 	
		);
int q2e(const float *pq,
		float *pe);
int e2r(const float *pe, float *pr);
int rvec2q(const float *prvec, float *pq);
int q2rvec(const float *pq, float *prvec);
int rvec2e(float *prvec, float *pe);
void liftPt(const _CamCoeff camCoeff, float *srcPt, float *dstPt);
void cam2imu(const float *pR, const float *pT, float *srcPt, float *dstPt);
void b2n(const float *pq, float *srcPt, float *dstPt);
void n2b(const float *pq, float *srcPt, float *dstPt);
void n2h(const float *pq, float *srcPt, float *dstPt);
void transpose(const float *pr, float *prt);
void h2n(const float *pq, float *srcPt, float *dstPt);
void b2h(const float *pq, float *srcPt, float *dstPt);
void h2b(const float *pq, float *srcPt, float *dstPt);


#endif
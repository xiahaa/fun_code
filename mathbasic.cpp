#include "mathbasic.h"
#include "opencv/cv.h"
#include <math.h>


int q2r(const float *pq,   //1x4 quaternions 
		float *pr	   //3x3 rotation matrix	
		)
{
	float q0, q1, q2, q3;
	q0 = pq[0]; q1 = pq[1]; q2 = pq[2]; q3 = pq[3];
	pr[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
	pr[1] = 2*(q1*q2 - q0*q3);
	pr[2] = 2*(q1*q3 + q0*q2);
	pr[3] = 2*(q1*q2 + q0*q3);
	pr[4] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
	pr[5] = 2*(q2*q3 - q0*q1);
	pr[6] = 2*(q1*q3 - q0*q2);
	pr[7] = 2*(q2*q3 + q0*q1);
	pr[8] = q0*q0 - q1*q1 - q2*q2 + q3*q3;

	return 0;
}

int r2e(const float *pr,   //3x3 rotation matrix	 
		float *pe	   //3x1 euler angles, 	
		)
{
	float roll, pitch, yaw;
	yaw = atan2(pr[3],pr[0]);
	pitch = asin(-pr[6]);
	roll = atan2(pr[7],pr[8]);
	
	pe[0] = yaw;
	pe[1] = pitch;
	pe[2] = roll;

	return 1;
}

int q2e(const float *pq,
		float *pe)
{
	//float pr[3][3];
	
	//q2r(pq,pr);
	//r2e(pr,pe);

	float q0, q1, q2, q3;
	q0 = pq[0]; q1 = pq[1]; q2 = pq[2]; q3 = pq[3];

	float pr[9];
	
	pr[0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
	pr[3] = 2*(q1*q2 + q0*q3);
	pr[6] = 2*(q1*q3 - q0*q2);
	pr[7] = 2*(q2*q3 + q0*q1);
	pr[8] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
	
	pe[0] = atan2(pr[3],pr[0]);
	pe[1] = asin(-pr[6]);
	pe[2] = atan2(pr[7],pr[8]);

	return 1;
}

int e2r(const float *pe, float *pr)
{
	float y = pe[0];
	float p = pe[1];
	float r = pe[2];

	pr[0] = cos(p)*cos(y);
	pr[1] = cos(p)*sin(y);
	pr[2] = -sin(p);
	pr[3] = sin(r)*sin(p)*cos(y) - cos(r)*sin(y);
	pr[4] = sin(r)*sin(p)*sin(y) + cos(r)*cos(y);
	pr[5] = cos(p)*sin(r);
	pr[6] = cos(r)*sin(p)*cos(y) + sin(r)*sin(y);
	pr[7] = cos(r)*sin(p)*sin(y) - sin(r)*cos(y);
	pr[8] = cos(p)*cos(r);


	/*float rrx[] = {
		1,		0,		 0,
		0,	cos(r),	sin(r),
		0, -sin(r),	cos(r),
	};
	float rry[] = {
	cos(p),		0, -sin(p),
		0,		1,		 0,
	sin(p),		0,	cos(p),
	};
	float rrz[] = {
	cos(y),	sin(y),		0,
   -sin(y),	cos(y),		0,
		0,		 0,		1,
	};
	
	CvMat hmat_Rx = cvMat(3, 3, CV_32FC1, rrx);
	CvMat hmat_Ry = cvMat(3, 3, CV_32FC1, rry);
	CvMat hmat_Rz = cvMat(3, 3, CV_32FC1, rrz);
	CvMat mat_rot = cvMat(3, 3, CV_32FC1, pr);

	cvMatMul(&hmat_Rx, &hmat_Ry, &mat_rot);
	cvMatMul(&mat_rot, &hmat_Rz, &mat_rot);*/

	return 0;
}

float norm3(const float *pvec)
{
	return sqrt(pvec[0] * pvec[0] + pvec[1] * pvec[1] + pvec[2] * pvec[2]);	
}

int rvec2q(const float *prvec, float *pq)
{
	float alpha_vec = norm3(prvec);
	
	pq[0] = cos(alpha_vec / 2.0);
	pq[1] = prvec[0] / alpha_vec * sin(alpha_vec / 2.0);
	pq[2] = prvec[1] / alpha_vec * sin(alpha_vec / 2.0);
	pq[3] = prvec[2] / alpha_vec * sin(alpha_vec / 2.0);

	return 1;
}

int q2rvec(const float *pq, float *prvec)
{
	float q0, q1, q2, q3;
	q0 = pq[0]; q1 = pq[1]; q2 = pq[2]; q3 = pq[3];
	float alpha = 2 * acos(q0);
	float norm_vec = norm3(pq+1);
	prvec[0] = alpha * q1 / norm_vec;
	prvec[1] = alpha * q2 / norm_vec;
	prvec[2] = alpha * q3 / norm_vec;

	return 1;
}

int rvec2e(float *prvec, float *pe)
{
	float pq[4];
	rvec2q(prvec,pq);
	q2e(pq,pe);
	return 1;
}

void liftPt(const _CamCoeff camCoeff, float *srcPt, float *dstPt)
{
	float cx = camCoeff.cx;
	float cy = camCoeff.cy;
	float ux = camCoeff.ux;
	float uy = camCoeff.uy;

	dstPt[0] = (srcPt[0] - cx) / ux;
	dstPt[1] = (srcPt[1] - cy) / uy;
	dstPt[2] = 1.f;
}

void cam2imu(const float *pR, const float *pT, float *srcPt, float *dstPt)
{
	dstPt[0] = pR[0] * srcPt[0] + pR[1] * srcPt[1] + pR[2] * srcPt[2] +	pT[0];
	dstPt[1] = pR[3] * srcPt[0] + pR[4] * srcPt[1] + pR[5] * srcPt[2] +	pT[1];
	dstPt[2] = pR[6] * srcPt[0] + pR[7] * srcPt[1] + pR[8] * srcPt[2] +	pT[2];
}

void b2n(const float *pq, float *srcPt, float *dstPt)
{
	float rnb[9];
	q2r(pq,rnb);
	dstPt[0] = rnb[0] * srcPt[0] + rnb[1] * srcPt[1] + rnb[2] * srcPt[2];
	dstPt[1] = rnb[3] * srcPt[0] + rnb[4] * srcPt[1] + rnb[5] * srcPt[2];
	dstPt[2] = rnb[6] * srcPt[0] + rnb[7] * srcPt[1] + rnb[8] * srcPt[2];
}

void n2b(const float *pq, float *srcPt, float *dstPt)
{
	float pe[3];
	q2e(pq,pe);
	float rbn[9];
	e2r(pe,rbn);
	dstPt[0] = rbn[0] * srcPt[0] + rbn[1] * srcPt[1] + rbn[2] * srcPt[2];
	dstPt[1] = rbn[3] * srcPt[0] + rbn[4] * srcPt[1] + rbn[5] * srcPt[2];
	dstPt[2] = rbn[6] * srcPt[0] + rbn[7] * srcPt[1] + rbn[8] * srcPt[2];
}

void n2h(const float *pq, float *srcPt, float *dstPt)
{
	float pe[3];
	q2e(pq,pe);
	float rhn[9];
	pe[2] = 0;
	pe[1] = 0;
	e2r(pe,rhn);
	dstPt[0] = rhn[0] * srcPt[0] + rhn[1] * srcPt[1] + rhn[2] * srcPt[2];
	dstPt[1] = rhn[3] * srcPt[0] + rhn[4] * srcPt[1] + rhn[5] * srcPt[2];
	dstPt[2] = rhn[6] * srcPt[0] + rhn[7] * srcPt[1] + rhn[8] * srcPt[2];
}

void transpose(const float *pr, float *prt)
{
	prt[0] = pr[0];
	prt[1] = pr[3];
	prt[2] = pr[6];
	prt[3] = pr[1];
	prt[4] = pr[4];
	prt[5] = pr[7];
	prt[6] = pr[2];
	prt[7] = pr[5];
	prt[8] = pr[8];
}


void h2n(const float *pq, float *srcPt, float *dstPt)
{
	float pe[3];
	q2e(pq,pe);
	float rhn[9];
	pe[2] = 0;
	pe[1] = 0;
	e2r(pe,rhn);
	float rnh[9];
	transpose(rhn, rnh);
	dstPt[0] = rnh[0] * srcPt[0] + rnh[1] * srcPt[1] + rnh[2] * srcPt[2];
	dstPt[1] = rnh[3] * srcPt[0] + rnh[4] * srcPt[1] + rnh[5] * srcPt[2];
	dstPt[2] = rnh[6] * srcPt[0] + rnh[7] * srcPt[1] + rnh[8] * srcPt[2];
}

void b2h(const float *pq, float *srcPt, float *dstPt)
{
	float tmpPt[3];
	b2n(pq,srcPt,tmpPt);
	n2h(pq,tmpPt,dstPt);
}

void h2b(const float *pq, float *srcPt, float *dstPt)
{
	float tmpPt[3];
	h2n(pq,srcPt,tmpPt);
	n2b(pq,tmpPt,dstPt);
}
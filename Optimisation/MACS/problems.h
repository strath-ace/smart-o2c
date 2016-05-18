// cec09.h
//
// C++ source codes of the test instances for CEC 2009 MOO competition
//
// If the source codes are not consistant with the report, please use the version in the report
//
// History
// 	  v1   Sept.8  2008
// 	  v1.1 Sept.22 2008: add R2_DTLZ2_M5 R3_DTLZ3_M5 WFG1_M5
//    v1.2 Oct.2  2008: fix the bugs in CF1-CF4, CF6-CF10, thank Qu Boyang for finding the bugs
//    v1.3 Oct.8  2008: fix a bug in R2_DTLZ2_M5, thank Santosh Tiwari
//    v1.4 Nov.26 2008: fix the bugs in CF4-CF5, thank Xueqiang Li for finding the bugs
//    v1.5 Dec.2  2008: fix the bugs in CF5 and CF7, thank Santosh Tiwari for finding the bugs
	
#ifndef AZ_CEC09_H
#define AZ_CEC09_H

namespace CEC09
{	
	void R2_DTLZ2_M5(double *x, double *f, const unsigned int nx, const unsigned int n_obj);
	void R3_DTLZ3_M5(double *x, double *f, const unsigned int nx, const unsigned int n_obj);
	void WFG1_M5( double *z, double *f, const unsigned int nx,  const unsigned int M);
}
#endif //AZ_CEC09_H

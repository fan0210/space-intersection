#ifndef SGM_H_
#define SGM_H_

#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
using namespace cv;

typedef float SGM_DP_TYPE;

class SGM
{
public:
	SGM();
	~SGM();

	//得到CT代价
	int Count1(unsigned int x);
	void CalulateCT(unsigned char *src, int width, int height, unsigned int *dst);
	void hamin(unsigned int *src1, unsigned *src2, int DMAX, int width, int height, float lamada, vector< vector<int> > &cpds);
	//得到SAD代价
	void SAD(uchar *src1, uchar *src2, int height, int width, int DMAX, int WindowSize, vector<vector<float> > &cpds);

	//动态规划相关
	int GetP2(unsigned char G);
	SGM_DP_TYPE min4(SGM_DP_TYPE d1, SGM_DP_TYPE d2, SGM_DP_TYPE d3, SGM_DP_TYPE d4);
	void SGM_DP_L(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_R(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_LU(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_RD(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_RU(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_LD(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_U(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_D(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp, unsigned char *img_data, vector<vector<SGM_DP_TYPE>>& new_cpds);
	void SGM_DP_avg2(vector< vector<SGM_DP_TYPE> >& cpds_1, vector<vector<SGM_DP_TYPE>>&cpds_2, vector<vector<SGM_DP_TYPE>>&cpds_12, int width, int height, int disp_N);

	void SGM_DP_min2(vector< vector<SGM_DP_TYPE> >& cpds_1, vector<vector<SGM_DP_TYPE>>&cpds_2, vector<vector<SGM_DP_TYPE>>&cpds_12, int width, int height, int disp_N);

	void SGM_DP_avg4(vector< vector<SGM_DP_TYPE> >& cpds_1, vector<vector<SGM_DP_TYPE>>&cpds_2, vector< vector<SGM_DP_TYPE> >& cpds_3, vector<vector<SGM_DP_TYPE>>&cpds_4, vector<vector<SGM_DP_TYPE>>&cpds_1234, int width, int height, int disp_N);
	void SGM_DP_avg5(vector< vector<SGM_DP_TYPE> >& cpds_1, vector<vector<SGM_DP_TYPE>>&cpds_2, vector< vector<SGM_DP_TYPE> >& cpds_3, vector<vector<SGM_DP_TYPE>>&cpds_4, vector<vector<SGM_DP_TYPE>>&cpds_5, vector<vector<SGM_DP_TYPE>>&cpds_12345, int width, int height, int disp_N);
	void SGM_DP_avg8(vector< vector<SGM_DP_TYPE> >& cpds_1, vector<vector<SGM_DP_TYPE>>&cpds_2, vector< vector<SGM_DP_TYPE> >& cpds_3, vector<vector<SGM_DP_TYPE>>&cpds_4, vector< vector<SGM_DP_TYPE> >& cpds_5, vector<vector<SGM_DP_TYPE>>&cpds_6, vector< vector<SGM_DP_TYPE> >& cpds_7, vector<vector<SGM_DP_TYPE>>&cpds_8, vector<vector<SGM_DP_TYPE>>&cpds_o, int width, int height, int disp_N);
	//得到视差图
	void GetDisp(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height, int disp_N, unsigned int *disp);

	//LRcheck
	template <typename T>
	void GetRightDisp(unsigned int* rightdisp, vector<vector<T >> &cost, int width, int height, int DMax)
	{


		for (int i = 0; i < height; i++)
			for (int j = 0; j < width; j++)

			{
				/*int destdisp = 0;
				int mincost = cost[i*width + j][0];
				int secmincost = cost[i*width + j][0];*/
				int maxcurrdisp = ((DMax - 1) - (width - 1 - j)) < 0 ? (DMax - 1) : (width - 1 - j);
				vector<float> rightcost;
				if (maxcurrdisp == (DMax - 1))
				{
					for (int d = 0; d < DMax; d++)

					{

						//cout << cost[i*width + j + d][d] << endl;


						T tmp = cost[i*width + j + d][d];
						rightcost.push_back(tmp);

					}

					vector< float >::iterator minist = std::min_element(std::begin(rightcost), std::end(rightcost));
					unsigned char min_index = minist - std::begin(rightcost);
					rightdisp[i*width + j] = min_index;
				}
				else
				{
					int destdisp = 0;
					int mincost = cost[i*width + j][0];
					int secmincost = cost[i*width + j][0];

					for (int k = 0; k < maxcurrdisp - 1; k++)

					{
						rightcost.push_back(cost[i*width + j + k][k]);
					}
					vector< float >::iterator minist1 = std::min_element(std::begin(rightcost), std::end(rightcost));
					unsigned char min_index = minist1 - std::begin(rightcost);
					rightdisp[i*width + j] = min_index;
				}
			}
	}
	template <typename T>
	void LRcheck(int height, int width, T *leftdisp, T *rightdisp, T *LRcheckresult, int DMax)
	{
		for (int row = 0; row < height; row++)
			for (int col = 0; col < width; col++)

			{
				//cout << leftdisp[0] << "   " << rightdisp[0];
				int leftpixelDisp = leftdisp[row*width + col];
				if (row*width + col - leftpixelDisp >= 0)//小于0也是误匹配点
				{
					int rightpixelDisp = rightdisp[row*width + col - leftpixelDisp];

					//	cout << leftdisp[0]<<"   "<< rightdisp[0];
					int diff = rightpixelDisp - leftpixelDisp;

					if (abs(diff) > 1)
					{
						LRcheckresult[row*width + col] = 0;//标志误匹配点
														   //	cout << row << "  " << col << endl;

					}
					else
					{
						LRcheckresult[row*width + col] = leftpixelDisp;
					}
				}
				else
					LRcheckresult[row*width + col] = 1;//标志误匹配点
			}

	}
	//亚像素插值
	void subPixelrefine(const Mat &LRcheckout, vector<vector< SGM_DP_TYPE>> spds, Mat &refineout);
	//显示
	void GenerateFalseMap(cv::Mat &src, cv::Mat &disp);
	void ShowcolorDisp(unsigned int *disp, int width, int height, char * name, char *savename, int DMax);
	void ShowDispgray(unsigned int *disp, int width, int height, char * name, char *savename, int DMax);

	//输入图像参数不同的立体匹配接口

	int GetDisprity_mat(cv::Mat& img11, cv::Mat &img22, int path, int DMax, bool debug_display, Mat &finaldisp);
	int Getdisprity_float(float*left, float*right, int path, int DMax, bool debug_display, int width, int height, float*lrcheckdisp, float*rawdisp);

	//其他工具
	void uintToMat_8UC1(unsigned int *input, cv::Mat& out, int height, int width);
	void uintTofloat(unsigned int *input, float* out, int height, int width);
private:
	static const	float P1;
	static const float P2_apha;
	static const float P2_gamma;
	static const float P2_min;
};

#endif
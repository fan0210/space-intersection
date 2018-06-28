#include"SGM.h"

 const	float SGM::P1 = 7;
 const float  SGM::P2_apha  = 0.25f;
 const float SGM::P2_gamma = 100;
 const float SGM::P2_min   = 17;
 SGM::SGM(){}
 SGM::~SGM(){}
//得到CT代价
int SGM::Count1(unsigned int x)
{
	int c = 0;
	while (x)
	{
		x &= x - 1;
		c++;
	}
	return c;
}
void SGM::CalulateCT(unsigned char *src, int width, int height, unsigned int *dst)
{
	int index = 0;
	memset(src, width*height, 0);
	for (int row = 2; row < height - 2; row++)
		for (int col = 2; col < width - 2; col++)
		{
			int res = 0;
			int midval = src[(row)*width + col];
			for (int m = -2; m <= 2; m++)
				for (int n = -2; n <= 2; n++)
				{
					if (m == 0 && n == 0) continue;
					res <<= 1;
					if (src[(row + m)*width + col + n] > midval)
					{
						res |= 1;
					}
				}
			dst[row*width + col] = res;
		}
}
void SGM::hamin(unsigned int *src1, unsigned *src2, int DMAX, int width, int height, float lamada, vector< vector<int> > &cpds)
{ 
	for (int row = 0; row < height; row++)
	{
		int ii = row;//没用到
		for (int col = 0; col < width; col++)
		{
			vector<int> cpd;
			unsigned int left = src1[row*width + col];
			for (int i = 0; i < DMAX; i++)
			{
				if (col >= i)
				{
					unsigned int right = src2[row*width + col - i];
					unsigned int xor_res = left^right;
					int cnt = Count1(xor_res);

					cpd.push_back(Count1(xor_res)>lamada ? lamada : Count1(xor_res));
					//of << cnt << " ";
				}
				else
				{
					cpd.push_back(12);
					//of << 12 << " ";
				}
			}
			cpds.push_back(cpd);
			//of << endl;
		}
	}
}
//得到SAD代价

void SGM::SAD(uchar *src1, uchar *src2, int height, int width, int DMAX, int WindowSize, vector<vector<float> > &cpds)
{
	for (int row = 0; row < height; row++)
		for (int col = 0; col < width; col++)
		{
			vector<float> cpd;
			if (row <= WindowSize / 2 || col <= WindowSize / 2 || row >= height - WindowSize / 2 || col >= width - WindowSize / 2)//边界填充
			{
				for (int d = 0; d < DMAX; d++)
					cpd.push_back(1000);
			}
			else
			{
				for (int d = 0; d < DMAX; d++)
				{

					if (d < col - WindowSize / 2)
					{
						float sum = 0;
						for (int i = -WindowSize / 2; i <= WindowSize / 2; i++)
							for (int j = -WindowSize / 2; j <= WindowSize / 2; j++)
							{
								float left = src1[(row + i)*width + col + j];
								float right = src2[(row + i)*width + col + j - d];

								sum += abs(left - right);
							}
						cpd.push_back(sum);

					}
					else
					{
						cpd.push_back(1000);

					}

				}
			}

			cpds.push_back(cpd);
		}
}
 //动态规划相关
int SGM::GetP2(unsigned char G)
{
	 int result;
     result = (int)(-P2_apha *G + P2_gamma);
    if (result < P2_min)
        result = P2_min;
	return result;
	
}
SGM_DP_TYPE SGM::min4(SGM_DP_TYPE d1,SGM_DP_TYPE d2,SGM_DP_TYPE d3,SGM_DP_TYPE d4)
{
	SGM_DP_TYPE min12 = d1 <= d2 ? d1 : d2;
	SGM_DP_TYPE min34 = d3 <= d4 ? d3 : d4;

	return min12 < min34 ? min12 : min34;
}
void SGM::SGM_DP_L(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N, unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)
{
	int invalid_data = 0;
	vector<SGM_DP_TYPE> Lr_Left(cpds[0]);

	SGM_DP_TYPE Lr1_Left_min       = invalid_data;
	int     Lr1_Left_min_index = 0;


	for(int row = 0; row < height; row ++)
	for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd = cpds[row*width + col];
		vector<SGM_DP_TYPE> Lr;
		if(col==0)
		{
			Lr = cpd;
		}else //update
		{
			unsigned char G ;
			if(col < 0)
					G = 0;
			else
				G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[row*width + col -1]));
			int P2 = GetP2(G);
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_Left[u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_Left[u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_Left[u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr1_Left_min +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr1_Left_min);
			}
		}
		new_cpds.push_back(Lr);

		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
		Lr_Left = Lr;
		Lr1_Left_min = *minist;
		Lr1_Left_min_index = min_index;
		if(disp)
			disp[row*width + col] =  min_index;
	}
}

void SGM::SGM_DP_R(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)
{

	//
	new_cpds = vector<vector<SGM_DP_TYPE>>(width*height,vector<SGM_DP_TYPE>(disp_N,0));
	int invalid_data = 0;
	vector<SGM_DP_TYPE> Lr_Right;

	SGM_DP_TYPE Lr_Right_min       = invalid_data;
	int     Lr_Right_min_index = 0;


	for(int row = 0; row < height; row ++)
	for(int col = width-1; col >= 0 ; col --)
	{
		vector<SGM_DP_TYPE> cpd = cpds[row*width + col];
		vector<SGM_DP_TYPE> Lr;
		if(col== width -1)
		{
			Lr = cpd;
		}else
		{
			int G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[row*width + col + 1]));

			int P2 = GetP2(G);
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_Right[u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_Right[u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_Right[u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr_Right_min +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr_Right_min);
			}
		}
		new_cpds[row*width + col] = (Lr);

		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
		Lr_Right = Lr;
		Lr_Right_min = *minist;
		Lr_Right_min_index = min_index;
		if(disp)
			disp[row*width + col] =  min_index;
	}
}

void SGM::SGM_DP_LU(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)
{
	int invalid_data1 = 0;
	vector< vector<SGM_DP_TYPE> > Lr_UPs(width, vector<SGM_DP_TYPE>(64,invalid_data1));     
	vector<SGM_DP_TYPE> Lr_UPs_min(width,invalid_data1);
	vector<int> Lr_UPs_min_index(width,0);

	vector<SGM_DP_TYPE> Lr_UL;
	SGM_DP_TYPE         Lr_UL_min;
	int             Lr_UL_min_index;

	////第一行，Lr1已经可以运算
	for(int row = 0; row < height; row ++)
	for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd = cpds[row*width + col];
		vector<SGM_DP_TYPE> Lr;
		if(row == 0 || col == 0)
		{
			Lr = cpd;
		}else
		{
			int G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[(row-1)*width + col - 1]));
			int P2 = GetP2(G);
			
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_UL[u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_UL[u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_UL[u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr_UL_min +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr_UL_min);
			}
		}
		new_cpds.push_back(Lr);
		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
		Lr_UL     = Lr_UPs[col];
		Lr_UL_min = Lr_UPs_min[col];
		Lr_UL_min_index = Lr_UPs_min_index[col];

		Lr_UPs[col] = Lr;
		Lr_UPs_min[col] = *minist;
		Lr_UPs_min_index[col] = min_index;
		if(disp)
			disp[row*width + col] =  min_index;
	}

}

void SGM::SGM_DP_RD(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)
{
	new_cpds = vector<vector<SGM_DP_TYPE>>(width*height,vector<SGM_DP_TYPE>(disp_N,0));
	int invalid_data1 = 0;
	vector< vector<SGM_DP_TYPE> > Lr_UPs(width, vector<SGM_DP_TYPE>(64,invalid_data1));     
	vector<SGM_DP_TYPE> Lr_UPs_min(width,invalid_data1);
	vector<int> Lr_UPs_min_index(width,0);

	vector<SGM_DP_TYPE> Lr_UL;
	SGM_DP_TYPE         Lr_UL_min;
	int             Lr_UL_min_index;

	////第一行，Lr1已经可以运算
	for(int row = height - 1; row >= 0    ; row --)
	for(int col = width -1        ; col >= 0 ; col --)
	{
		vector<SGM_DP_TYPE> cpd = cpds[row*width + col];
		vector<SGM_DP_TYPE> Lr;
		if(row == height - 1 || col == width -1)
		{
			Lr = cpd;
		}else
		{
			int G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[(row+1)*width + col + 1]));
			int P2 = GetP2(G);
			
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_UL[u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_UL[u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_UL[u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr_UL_min +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr_UL_min);
			}
		}
		new_cpds[row*width + col] = (Lr);
		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
		Lr_UL     = Lr_UPs[col];
		Lr_UL_min = Lr_UPs_min[col];
		Lr_UL_min_index = Lr_UPs_min_index[col];

		Lr_UPs[col] = Lr;
		Lr_UPs_min[col] = *minist;
		Lr_UPs_min_index[col] = min_index;
		if(disp)
			disp[row*width + col] =  min_index;
	}

}

void SGM::SGM_DP_RU(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)
{
	
	int invalid_data1 = 0;
	vector< vector<SGM_DP_TYPE> > Lr_UPs(width, vector<SGM_DP_TYPE>(64,invalid_data1));     
	vector<SGM_DP_TYPE> Lr_UPs_min(width,invalid_data1);
	vector<int> Lr_UPs_min_index(width,0);

	vector<SGM_DP_TYPE> Lr_UR;
	SGM_DP_TYPE         Lr_UR_min;
	int             Lr_UR_min_index;

	////第一行，Lr1已经可以运算
	for(int row = 0; row < height; row ++)
	for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd = cpds[row*width + col]; //获取当前代价
		vector<SGM_DP_TYPE> Lr;
		if(col < width-1)
		{
			//对应右上的位置
			Lr_UR     = Lr_UPs[col+1];
			Lr_UR_min = Lr_UPs_min[col+1];
			Lr_UR_min_index = Lr_UPs_min_index[col+1];
		}
		if(row == 0 || col == (width-1))
		{
			Lr = cpd;
		}else
		{
			int G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[(row+1)*width + col + 1]));
			int P2 = GetP2(G);
			
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_UR[u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_UR[u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_UR[u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr_UR_min +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr_UR_min);
			}
		}
		new_cpds.push_back(Lr);
		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
	
		Lr_UPs[col] = Lr;
		Lr_UPs_min[col] = *minist;
		Lr_UPs_min_index[col] = min_index;
		if(disp)
			disp[row*width + col] =  min_index;
	}

}

void SGM::SGM_DP_LD(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)

{
	new_cpds = vector<vector<SGM_DP_TYPE>>(width*height,vector<SGM_DP_TYPE>(disp_N,0));
	int invalid_data1 = 0;
	vector< vector<SGM_DP_TYPE> > Lr_UPs(width, vector<SGM_DP_TYPE>(64,invalid_data1));     
	vector<SGM_DP_TYPE> Lr_UPs_min(width,invalid_data1);
	vector<int> Lr_UPs_min_index(width,0);

	vector<SGM_DP_TYPE> Lr_UR;
	SGM_DP_TYPE         Lr_UR_min;
	int             Lr_UR_min_index;

	////第一行，Lr1已经可以运算
	for(int row = height - 1; row >= 0; row --)
	for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd = cpds[row*width + col];
		vector<SGM_DP_TYPE> Lr;
		if(col > 0)
		{
			Lr_UR     = Lr_UPs[col-1];
			Lr_UR_min = Lr_UPs_min[col-1];
			Lr_UR_min_index = Lr_UPs_min_index[col-1];
		}
		if(row == height - 1 || col == 0)
		{
			Lr = cpd;
		}else
		{
			int G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[(row+1)*width + col - 1]));
			int P2 = GetP2(G);
			
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_UR[u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_UR[u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_UR[u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr_UR_min +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr_UR_min);
			}
		}
		new_cpds[row*width + col] = (Lr);
		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
		Lr_UPs[col] = Lr;
		Lr_UPs_min[col] = *minist;
		Lr_UPs_min_index[col] = min_index;
		if(disp)
			disp[row*width + col] =  min_index;
	}

}

void SGM::SGM_DP_U(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)
{
	int invalid_data1 = 0;
	vector< vector<SGM_DP_TYPE> > Lr_Ups(width, vector<SGM_DP_TYPE>(disp_N,invalid_data1));      vector<SGM_DP_TYPE> Lr_Ups_mins(width,invalid_data1); vector<int> Lr_Ups_mins_index(width,0);

	////第一行，Lr1已经可以运算
	//for(int col = 1; col <)
	for(int row = 0; row < height; row ++)
	for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE>& cpd = cpds[row*width + col];
		vector<SGM_DP_TYPE>  Lr;
		if(row == 0)
		{
			Lr = cpd;
		}else
		{
			int G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[(row-1)*width + col]));
			int P2 = GetP2(G);
			
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_Ups[col][u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_Ups[col][u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_Ups[col][u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr_Ups_mins[col] +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr_Ups_mins[col]);
			}
		}
		new_cpds.push_back(Lr);
		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
		Lr_Ups[col] = Lr;
		Lr_Ups_mins[col] = *minist;
		Lr_Ups_mins_index[col] = min_index;
		if(disp)
			disp[row*width + col] =  min_index;
		
	}
}
void SGM::SGM_DP_D(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp, unsigned char *img_data,vector<vector<SGM_DP_TYPE>>& new_cpds)
{

	new_cpds = vector<vector<SGM_DP_TYPE>>(width*height,vector<SGM_DP_TYPE>(disp_N,0));
	int invalid_data1 = 0;
	vector< vector<SGM_DP_TYPE> > Lr_Ups(width, vector<SGM_DP_TYPE>(64,invalid_data1));      vector<SGM_DP_TYPE> Lr_Ups_mins(width,invalid_data1); vector<int> Lr_Ups_mins_index(width,0);

	////第一行，Lr1已经可以运算
	//for(int col = 1; col <)
	for(int row = height - 1; row >= 0; row --)
	for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE>& cpd = cpds[row*width + col];
		vector<SGM_DP_TYPE>  Lr;
		if(row == height - 1)
		{
			Lr = cpd;
		}else
		{
			int G = abs((int)((unsigned int)img_data[row*width + col]) - (int)((unsigned int)img_data[(row+1)*width + col]));
			int P2 = GetP2(G);
			
			for(int u = 0; u < disp_N; u++)
			{
				SGM_DP_TYPE d_elem        = Lr_Ups[col][u];
				SGM_DP_TYPE d_s1_elem     = u== 0        ? UCHAR_MAX: Lr_Ups[col][u-1] + P1;
				SGM_DP_TYPE d_p1_elem     = u== disp_N-1 ? UCHAR_MAX: Lr_Ups[col][u+1] + P1;
				SGM_DP_TYPE min_elem      = Lr_Ups_mins[col] +  P2;
				SGM_DP_TYPE minist = min4(d_elem,d_s1_elem,d_p1_elem,min_elem);
				Lr.push_back(minist  + cpd[u] - Lr_Ups_mins[col]);
			}
		}
		new_cpds[row*width + col] = (Lr);
		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(Lr), std::end(Lr));  
		int min_index = minist - std::begin(Lr);
		Lr_Ups[col] = Lr;
		Lr_Ups_mins[col] = *minist;
		Lr_Ups_mins_index[col] = min_index;
		if(disp)
		disp[row*width + col] =  min_index;
	}
}
void SGM::SGM_DP_avg2(vector< vector<SGM_DP_TYPE> >& cpds_1,vector<vector<SGM_DP_TYPE>>&cpds_2,vector<vector<SGM_DP_TYPE>>&cpds_12, int width, int height,int disp_N)
{
	for(int row = 0; row < height; row ++)
		for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd_1 = cpds_1[row*width + col];
		vector<SGM_DP_TYPE> cpd_2 = cpds_2[row*width + col];
		vector<SGM_DP_TYPE> cpd_12;
		for(int u = 0; u < disp_N; u++)
		{
			SGM_DP_TYPE cpd_1_u = cpd_1[u];
			SGM_DP_TYPE cpd_2_u = cpd_2[u];
			cpd_12.push_back((cpd_1_u+cpd_2_u)/2);
		}
		cpds_12.push_back(cpd_12);
	}
}

void SGM::SGM_DP_min2(vector< vector<SGM_DP_TYPE> >& cpds_1,vector<vector<SGM_DP_TYPE>>&cpds_2,vector<vector<SGM_DP_TYPE>>&cpds_12, int width, int height,int disp_N)
{
	//vector< vector<SGM_DP_TYPE> > cpds_LR;
	for(int row = 0; row < height; row ++)
		for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd_1 = cpds_1[row*width + col];
		vector<SGM_DP_TYPE> cpd_2 = cpds_2[row*width + col];
		vector<SGM_DP_TYPE> cpd_12;
		for(int u = 0; u < disp_N; u++)
		{
			SGM_DP_TYPE cpd_1_u = cpd_1[u];
			SGM_DP_TYPE cpd_2_u = cpd_2[u];
			cpd_12.push_back(cpd_1_u<=cpd_2_u?cpd_1_u:cpd_2_u);
			//cpd_LR.push_back((cpd_L_u+cpd_R_u)/2);
		}
		cpds_12.push_back(cpd_12);
	}
}

void SGM::SGM_DP_avg4(vector< vector<SGM_DP_TYPE> >& cpds_1,vector<vector<SGM_DP_TYPE>>&cpds_2,vector< vector<SGM_DP_TYPE> >& cpds_3,vector<vector<SGM_DP_TYPE>>&cpds_4,vector<vector<SGM_DP_TYPE>>&cpds_1234, int width, int height,int disp_N)
{
	//vector< vector<SGM_DP_TYPE> > cpds_1R;
	for(int row = 0; row < height; row ++)
		for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd_L = cpds_1[row*width + col];
		vector<SGM_DP_TYPE> cpd_R = cpds_2[row*width + col];
		vector<SGM_DP_TYPE> cpd_U = cpds_3[row*width + col];
		vector<SGM_DP_TYPE> cpd_D = cpds_4[row*width + col];
		vector<SGM_DP_TYPE> cpd_UD;
		for(int u = 0; u < disp_N; u++)
		{
			SGM_DP_TYPE cpd_L_u = cpd_L[u];
			SGM_DP_TYPE cpd_R_u = cpd_R[u];
			SGM_DP_TYPE cpd_U_u = cpd_U[u];
			SGM_DP_TYPE cpd_D_u = cpd_D[u];
			cpd_UD.push_back((cpd_L_u*2+ cpd_R_u*2+ cpd_U_u+cpd_D_u)/6);
			//cpd_LR.push_back((cpd_L_u+cpd_R_u)/2);
		}
		cpds_1234.push_back(cpd_UD);
	}
}
void SGM:: SGM_DP_avg5(vector< vector<SGM_DP_TYPE> >& cpds_1,vector<vector<SGM_DP_TYPE>>&cpds_2,vector< vector<SGM_DP_TYPE> >& cpds_3,vector<vector<SGM_DP_TYPE>>&cpds_4,vector<vector<SGM_DP_TYPE>>&cpds_5,vector<vector<SGM_DP_TYPE>>&cpds_12345, int width, int height,int disp_N)
{
	//vector< vector<SGM_DP_TYPE> > cpds_1R;
	for(int row = 0; row < height; row ++)
		for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd_L = cpds_1[row*width + col];
		vector<SGM_DP_TYPE> cpd_R = cpds_2[row*width + col];
		vector<SGM_DP_TYPE> cpd_U = cpds_3[row*width + col];
		vector<SGM_DP_TYPE> cpd_D = cpds_4[row*width + col];
		vector<SGM_DP_TYPE> cpd_5 = cpds_5[row*width + col];
		vector<SGM_DP_TYPE> cpd_UD;
		for(int u = 0; u < disp_N; u++)
		{
			SGM_DP_TYPE cpd_L_u = cpd_L[u];
			SGM_DP_TYPE cpd_R_u = cpd_R[u];
			SGM_DP_TYPE cpd_U_u = cpd_U[u];
			SGM_DP_TYPE cpd_D_u = cpd_D[u];
			SGM_DP_TYPE cpd_5_u = cpd_5[u];
			cpd_UD.push_back((cpd_L_u+ cpd_R_u+ cpd_U_u+cpd_D_u+cpd_5_u)/5);
			//cpd_LR.push_back((cpd_L_u+cpd_R_u)/2);
		}
		cpds_12345.push_back(cpd_UD);
	}
}

void SGM::SGM_DP_avg8(vector< vector<SGM_DP_TYPE> >& cpds_1,vector<vector<SGM_DP_TYPE>>&cpds_2,vector< vector<SGM_DP_TYPE> >& cpds_3,vector<vector<SGM_DP_TYPE>>&cpds_4,vector< vector<SGM_DP_TYPE> >& cpds_5,vector<vector<SGM_DP_TYPE>>&cpds_6,vector< vector<SGM_DP_TYPE> >& cpds_7,vector<vector<SGM_DP_TYPE>>&cpds_8,vector<vector<SGM_DP_TYPE>>&cpds_o, int width, int height,int disp_N)
{
	//vector< vector<SGM_DP_TYPE> > cpds_LR;
	for(int row = 0; row < height; row ++)
		for(int col = 0; col < width ; col ++)
	{
		vector<SGM_DP_TYPE> cpd_1 = cpds_1[row*width + col];
		vector<SGM_DP_TYPE> cpd_2 = cpds_2[row*width + col];
		vector<SGM_DP_TYPE> cpd_3 = cpds_3[row*width + col];
		vector<SGM_DP_TYPE> cpd_4 = cpds_4[row*width + col];
		vector<SGM_DP_TYPE> cpd_5 = cpds_5[row*width + col];
		vector<SGM_DP_TYPE> cpd_6 = cpds_6[row*width + col];
		vector<SGM_DP_TYPE> cpd_7 = cpds_7[row*width + col];
		vector<SGM_DP_TYPE> cpd_8 = cpds_8[row*width + col];
		vector<SGM_DP_TYPE> cpd_o;
		for(int u = 0; u < disp_N; u++)
		{
			SGM_DP_TYPE cpd_1_u = cpd_1[u];
			SGM_DP_TYPE cpd_2_u = cpd_2[u];
			SGM_DP_TYPE cpd_3_u = cpd_3[u];
			SGM_DP_TYPE cpd_4_u = cpd_4[u];
			SGM_DP_TYPE cpd_5_u = cpd_1[u];
			SGM_DP_TYPE cpd_6_u = cpd_2[u];
			SGM_DP_TYPE cpd_7_u = cpd_3[u];
			SGM_DP_TYPE cpd_8_u = cpd_4[u];
			cpd_o.push_back((cpd_1_u+ cpd_2_u+ cpd_3_u+cpd_4_u+cpd_5_u+ cpd_6_u+ cpd_7_u+cpd_8_u)/8);
			//cpd_LR.push_back((cpd_L_u+cpd_R_u)/2);
		}
		cpds_o.push_back(cpd_o);
	}
}


 //得到视差图
void SGM::GetDisp(vector< vector<SGM_DP_TYPE> >& cpds, int width, int height,int disp_N,unsigned int *disp)
{

	for(int row = 0; row < height; row ++)
	for(int col = 0; col <width  ; col ++)
	{
		vector<SGM_DP_TYPE> cpd_LR = cpds[row*width + col];;
		vector<SGM_DP_TYPE>::iterator minist = std::min_element(std::begin(cpd_LR), std::end(cpd_LR));  
		int min_index = minist - std::begin(cpd_LR);
		if(disp)
			disp[row*width + col] =  min_index;
	}

}
//亚像素插值
void SGM::subPixelrefine(const Mat &LRcheckout,vector<vector< SGM_DP_TYPE>> spds,Mat &refineout)
{
	int height=LRcheckout.rows;
	int width=LRcheckout.cols;
	refineout=cv::Mat::zeros(LRcheckout.size(),CV_32FC1);
            for(int row=0;row<height;row++)
				for(int col=0;col<width;col++)
				{
				     if(LRcheckout.at<uchar>(row,col)==0) continue;
					 else
					 {          
					         int d= LRcheckout.at<uchar>(row,col);
							SGM_DP_TYPE Spd=spds[row*width+col][d];
					        SGM_DP_TYPE Spd_last=spds[row*width+col][d-1];
							SGM_DP_TYPE Spd_next=spds[row*width+col][d+1];
							refineout.at<float>(row,col)=d+(Spd_last-Spd_next)/(2*Spd_next+2*Spd_last-4*Spd);
							
							//cout<<"LRcheckout"<<(int)LRcheckout.at<uchar>(row,col)<<endl;
						//	cout<<"refineout"<<refineout.at<float>(row,col)<<endl;
					 }				
				}
}

//显示
void SGM::GenerateFalseMap(cv::Mat &src, cv::Mat &disp)
{
    // color map
    float max_val = 255.0f;
    float map[8][4] = {{0,0,0,114},{0,0,1,185},{1,0,0,114},{1,0,1,174},
                       {0,1,0,114},{0,1,1,185},{1,1,0,114},{1,1,1,0}};
    float sum = 0;
    for (int i=0; i<8; i++)
      sum += map[i][3];

    float weights[8]; // relative   weights
    float cumsum[8];  // cumulative weights
    cumsum[0] = 0;
    for (int i=0; i<7; i++) {
      weights[i]  = sum/map[i][3];
      cumsum[i+1] = cumsum[i] + map[i][3]/sum;
    }

    int height_ = src.rows;
    int width_ = src.cols;
    // for all pixels do
    for (int v=0; v<height_; v++) {
      for (int u=0; u<width_; u++) {

        // get normalized value
        float val = std::min(std::max(src.data[v*width_ + u]/max_val,0.0f),1.0f);

        // find bin
        int i;
        for (i=0; i<7; i++)
          if (val<cumsum[i+1])
            break;

        // compute red/green/blue values
        float   w = 1.0-(val-cumsum[i])*weights[i];
        uchar r = (uchar)((w*map[i][0]+(1.0-w)*map[i+1][0]) * 255.0);
        uchar g = (uchar)((w*map[i][1]+(1.0-w)*map[i+1][1]) * 255.0);
        uchar b = (uchar)((w*map[i][2]+(1.0-w)*map[i+1][2]) * 255.0);
		//rgb内存连续存放
        disp.data[v*width_*3 + 3*u + 0] = b;
        disp.data[v*width_*3 + 3*u + 1] = g;
        disp.data[v*width_*3 + 3*u + 2] = r;
      }
    }
}
void SGM::ShowcolorDisp(unsigned int *disp, int width, int height, char * name, char *savename, int DMax)
{
	//cv::Mat a;//创建矩阵头

	//	a = cv::imread("f:\\psb.jpg");//读入图像
	IplImage * disp_orign = cvCreateImage(cvSize(width, height), 8, 1);//为IplImage型0-255视察图输出申请内存
	for (int row = 0; row < height; row++)
		for (int col = 0; col < width; col++)
		{
			float gray = disp[row*width + col] * float(256 / DMax);
			disp_orign->imageData[row*width + col] = gray;
		}
	///IplImage类型转Mat
	cv::Mat save = cv::Mat::zeros(disp_orign->height, disp_orign->width, CV_8UC1);
//	cv::cvarrToMat(disp_orign).copyTo(save);//IplImage类型（disp_orign）转Mat类型再复制到（save）里面
											//转化为彩色输出
	equalizeHist(cvarrToMat(disp_orign), save);//直方图均衡化
	cv::Mat colorresult = cv::Mat::zeros(height, width, CV_8UC3);//为彩色输出申请内存
	GenerateFalseMap(save, colorresult);//调用灰度转彩色函数
										//显示并保存
	imshow(name, colorresult);
	imwrite(savename, colorresult);

	//---------------------------------------------------------
	cvWaitKey(10);
	cvReleaseImage(&disp_orign);
}
void SGM::ShowDispgray(unsigned int *disp, int width, int height, char * name, char *savename, int DMax)
{
	//cv::Mat a;//创建矩阵头

	//	a = cv::imread("f:\\psb.jpg");//读入图像
	IplImage * disp_orign = cvCreateImage(cvSize(width, height), 8, 1);//为IplImage型0-255视察图输出申请内存
	for (int row = 0; row < height; row++)
		for (int col = 0; col < width; col++)
		{
			float gray = disp[row*width + col] * float(256 / DMax);
			disp_orign->imageData[row*width + col] = gray;
		}
	//Mat dst;
	//equalizeHist( Mat(disp_orign), dst );
	cvShowImage(name, disp_orign);
	cv::Mat save = cv::Mat::zeros(disp_orign->height, disp_orign->width, CV_8UC1);
	cv::cvarrToMat(disp_orign).copyTo(save);//IplImage类型（disp_orign）转Mat类型再复制到（save）里面
	imwrite(savename, save);
	//---------------------------------------------------------
	cvWaitKey(10);
	cvReleaseImage(&disp_orign);
}
//其他工具
void SGM::uintToMat_8UC1(unsigned int *input,cv::Mat& out ,int height,int width)
{ 
	for(int row=0;row<height;row++)
		for(int col=0;col<width;col++)
		{
		out.at<uchar>(row,col)=*(input+width*row+col);
		
		}


}
void SGM::uintTofloat(unsigned int *input,float *out ,int height,int width)
{
	for(int row=0;row<height;row++)
		for(int col=0;col<width;col++)
		{
		*(out+width*row+col)=*(input+width*row+col);
		
		}

}
//输入图像参数不同的立体匹配接口
 // 左图，右图，路径（1.2.4），是否显示得到的视差图，lrcheck后的视差图，未lrcheck后的视差图
int SGM::GetDisprity_mat(cv::Mat& img11,cv::Mat &img22,int path,int DMax,bool debug_display, Mat &finaldisp)

{
		             Mat LRcheckdisp,subrefineout;
					 Mat rawdisp;//先不用做输出，是以前结口的一部分
					  //对齐图像
					Mat img1_c_s, img2_c_s;
					Mat img1, img2;
					cv::resize(img11, img1_c_s, cv::Size(img11.cols, img11.rows));
					cv::resize(img22, img2_c_s, cv::Size(img22.cols, img22.rows));
				//	/////////////////
				
				switch (img1_c_s.cols%4)
				{
					
					case 1:
						{
							
							 img1 = Mat::zeros(img1_c_s.rows, img1_c_s.cols + 3, CV_8UC1);
							  img2 = Mat::zeros(img2_c_s.rows, img2_c_s.cols + 3, CV_8UC1);
						  for (int row = 0; row < img2_c_s.rows; row++)
							for (int col = 0; col < img2_c_s.cols - 3; col++)
							{
								img1.at<uchar>(row, col) = img1_c_s.data[img1_c_s.step*row + col];
								img2.at<uchar>(row, col) = img2_c_s.data[img1_c_s.step*row + col];
								
							}
							break;
						}

					case 2:
						{
							 img1 = Mat::zeros(img1_c_s.rows, img1_c_s.cols + 2, CV_8UC1);
							 img2 = Mat::zeros(img2_c_s.rows, img2_c_s.cols + 2, CV_8UC1);
						  for (int row = 0; row < img2_c_s.rows; row++)
							for (int col = 0; col < img2_c_s.cols - 2; col++)
							{
								img1.at<uchar>(row, col) = img1_c_s.data[img1_c_s.step*row + col];
								img2.at<uchar>(row, col) = img2_c_s.data[img1_c_s.step*row + col];
								
							}
							break;
						}
					case 3:
						{
							 img1 = Mat::zeros(img1_c_s.rows, img1_c_s.cols + 1, CV_8UC1);
							 img2 = Mat::zeros(img2_c_s.rows, img2_c_s.cols + 1, CV_8UC1);
						  for (int row = 0; row < img2_c_s.rows; row++)
							for (int col = 0; col < img2_c_s.cols - 1; col++)
							{
								img1.at<uchar>(row, col) = img1_c_s.data[img1_c_s.step*row + col];
								img2.at<uchar>(row, col) = img2_c_s.data[img1_c_s.step*row + col];
								
							}
							break;
						}
					default:
						{
						
							 img1 = img1_c_s;
							 img2=img2_c_s;
						
								break;
						
						}
				}
								

		if ((img1.rows&& img2.cols) == 0)
		{
			cout << " picture  is  not read in " << endl;
			system("pause");
			return -1;
		}
		
		int width = img1.cols;
		int height = img1.rows;

		
	

		unsigned int *disp = new unsigned int[width*height];
		
		uchar *src1 = img1.data;
		uchar *src2 = img2.data;
		//CT代价
			unsigned int *ct1 = new unsigned int[width*height];
			unsigned int *ct2 = new unsigned int[width*height];
			vector< vector<int> > cpds;
		
			CalulateCT(src1, width, height, ct1);
			CalulateCT(src2, width, height, ct2);
		
			hamin(ct1, ct2, DMax, width, height, 100, cpds);
		//SAD代价
		//vector< vector<float> > cpds;
		//SAD(src1,src2, height, width, DMax,3,cpds);		
		



		vector<vector<float>> costs_new;

		for (int i = 0; i < cpds.size(); i++)
		{
			vector<SGM_DP_TYPE> cost_new;
			for (int u = 0; u < DMax; u++)
			{
			
			cost_new.push_back(cpds[i][u]);
			}
			costs_new.push_back(cost_new);
		}
		
	
		
      switch (path)
	  {
		  	  

	  case 1:
		  {
			  vector<vector<SGM_DP_TYPE>> costs_L;
		  SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);															
		  GetDisp(costs_L, width, height, DMax, disp);
		   rawdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(disp, rawdisp , height, width);
		//先得到右视差图
			unsigned int *rightdisp=new unsigned int [width*height];
			 GetRightDisp(rightdisp,costs_L,  width,  height,DMax);
			 
			
		//再LRcheck
			unsigned int *LRcheckresult=new unsigned int [width*height];
			LRcheck(height,  width, disp, rightdisp, LRcheckresult, DMax);
			LRcheckdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(LRcheckresult, LRcheckdisp , height, width);
	
		  //亚像素插值
			
	        subPixelrefine(LRcheckdisp,costs_L,subrefineout);
		//中值滤波

			medianBlur(subrefineout,finaldisp,3);

					  if(debug_display)
					  {
						cv::Mat colorresult = cv::Mat::zeros(height, width, CV_8UC3);//为彩色输出申请内存
						cv::Mat tmp=LRcheckdisp.clone();
						medianBlur(LRcheckdisp,tmp,3);
						for(int row=0;row<height;row++)
							for(int col=0;col<width;col++)
							{
							tmp.at<uchar>(row,col)*=(float)256/DMax;//把视差值线性映射到0-255，只为显示。
							}
						  GenerateFalseMap(tmp,colorresult);
						  imshow(" ",colorresult);
						  imwrite("colordisp.png",colorresult);
					//  ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
					 //  ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
					   cvWaitKey(10);
					  }
		    delete  []rightdisp;
		    delete []  LRcheckresult;
		    delete []disp;
			break;  
		}
	  case 2:
		  {
			
			  vector<vector<SGM_DP_TYPE>> costs_L;
			  vector<vector<SGM_DP_TYPE>> costs_U;
			   vector<vector<SGM_DP_TYPE>> costs_LU;
		  SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);
		  SGM_DP_U(costs_new, width, height, DMax, disp, src1, costs_U); 
		  SGM_DP_avg2(costs_L, costs_U,  costs_LU, width, height, DMax);
		  GetDisp(costs_LU, width, height, DMax, disp);
		   rawdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(disp, rawdisp , height, width);
		//先得到右视差图
			unsigned int *rightdisp=new unsigned int [width*height];
			 GetRightDisp(rightdisp,costs_LU,  width,  height,DMax);
			 
			
		//再LRcheck
			unsigned int *LRcheckresult=new unsigned int [width*height];
			LRcheck(height,  width, disp, rightdisp, LRcheckresult, DMax);
			LRcheckdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(LRcheckresult, LRcheckdisp , height, width);
	
		  //亚像素插值
			
	        subPixelrefine(LRcheckdisp,costs_LU,subrefineout);
		//中值滤波

			medianBlur(subrefineout,finaldisp,3);

					  if(debug_display)
					  {
						cv::Mat colorresult = cv::Mat::zeros(height, width, CV_8UC3);//为彩色输出申请内存
						cv::Mat tmp=LRcheckdisp.clone();
						medianBlur(LRcheckdisp,tmp,3);
						for(int row=0;row<height;row++)
							for(int col=0;col<width;col++)
							{
							tmp.at<uchar>(row,col)*=(float)256/DMax;//把视差值线性映射到0-255，只为显示。
							}
						  GenerateFalseMap(tmp,colorresult);
						  imshow(" ",colorresult);
						  imwrite("colordisp.png",colorresult);
					//  ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
					 //  ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
					   cvWaitKey(10);
					  }
		
		      delete  []rightdisp;
		      delete []  LRcheckresult;
		      delete []disp;
			  break;
			  
		  }
		 case 4:
		  {
			
			  vector<vector<SGM_DP_TYPE>> costs_L;
			  vector<vector<SGM_DP_TYPE>> costs_U;
			  vector<vector<SGM_DP_TYPE>> costs_LU;
			  vector<vector<SGM_DP_TYPE>> costs_RU;
			   vector<vector<SGM_DP_TYPE>> AvgPath4;
		  SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);
		  SGM_DP_U(costs_new, width, height, DMax, disp, src1, costs_U); 
		  SGM_DP_LU(costs_new, width, height, DMax, disp, src1, costs_LU);
		  SGM_DP_RU(costs_new, width, height, DMax, disp, src1, costs_RU);
		  SGM_DP_avg4(costs_L, costs_U, costs_LU,costs_RU,AvgPath4, width, height, DMax);
		  GetDisp(AvgPath4, width, height, DMax, disp);
		   rawdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(disp, rawdisp , height, width);
		//先得到右视差图
			unsigned int *rightdisp=new unsigned int [width*height];
			 GetRightDisp(rightdisp,AvgPath4,  width,  height,DMax);
			 
			
		//再LRcheck
			unsigned int *LRcheckresult=new unsigned int [width*height];
			LRcheck(height,  width, disp, rightdisp, LRcheckresult, DMax);
			LRcheckdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(LRcheckresult, LRcheckdisp , height, width);
		//亚像素插值
			
	        subPixelrefine(LRcheckdisp,AvgPath4,subrefineout);
		//中值滤波

			medianBlur(subrefineout,finaldisp,3);

					  if(debug_display)
					  {
						cv::Mat colorresult = cv::Mat::zeros(height, width, CV_8UC3);//为彩色输出申请内存
						cv::Mat tmp=LRcheckdisp.clone();
						medianBlur(LRcheckdisp,tmp,3);
						for(int row=0;row<height;row++)
							for(int col=0;col<width;col++)
							{
							tmp.at<uchar>(row,col)*=(float)256/DMax;//把视差值线性映射到0-255，只为显示。
							}
						  GenerateFalseMap(tmp,colorresult);
						  imshow(" ",colorresult);
						  imwrite("colordisp.png",colorresult);
					//  ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
					 //  ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
					   cvWaitKey(10);
					  }
		      delete  []rightdisp;
		      delete []  LRcheckresult;
		      delete []disp;
			  break;
		  }
		   case 5:
		  {
			
			  vector<vector<SGM_DP_TYPE>> costs_L;
			  vector<vector<SGM_DP_TYPE>> costs_U;
			  vector<vector<SGM_DP_TYPE>> costs_LU;
			  vector<vector<SGM_DP_TYPE>> costs_RU;
			  vector<vector<SGM_DP_TYPE>> costs_R;
			   vector<vector<SGM_DP_TYPE>> AvgPath5;
		  SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);
		  SGM_DP_U(costs_new, width, height, DMax, disp, src1, costs_U); 
		  SGM_DP_LU(costs_new, width, height, DMax, disp, src1, costs_LU);
		  SGM_DP_RU(costs_new, width, height, DMax, disp, src1, costs_RU);
		   SGM_DP_R(costs_new, width, height, DMax, disp, src1, costs_R);
		  SGM_DP_avg5(costs_L, costs_U, costs_LU,costs_RU,costs_R,AvgPath5, width, height, DMax);
		  GetDisp(AvgPath5, width, height, DMax, disp);
		   rawdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(disp, rawdisp , height, width);
		//先得到右视差图
			unsigned int *rightdisp=new unsigned int [width*height];
			 GetRightDisp(rightdisp,AvgPath5,  width,  height,DMax);
			 
			
		//再LRcheck
			unsigned int *LRcheckresult=new unsigned int [width*height];
			LRcheck(height,  width, disp, rightdisp, LRcheckresult, DMax);
			LRcheckdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(LRcheckresult, LRcheckdisp , height, width);
			 //亚像素插值
			
	        subPixelrefine(LRcheckdisp,AvgPath5,subrefineout);
		//中值滤波

			medianBlur(subrefineout,finaldisp,3);

					  if(debug_display)
					  {
						cv::Mat colorresult = cv::Mat::zeros(height, width, CV_8UC3);//为彩色输出申请内存
						cv::Mat tmp=LRcheckdisp.clone();
						medianBlur(LRcheckdisp,tmp,3);
						for(int row=0;row<height;row++)
							for(int col=0;col<width;col++)
							{
							tmp.at<uchar>(row,col)*=(float)256/DMax;//把视差值线性映射到0-255，只为显示。
							}
						  GenerateFalseMap(tmp,colorresult);
						  imshow(" ",colorresult);
						  imwrite("colordisp.png",colorresult);
					//  ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
					 //  ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
					   cvWaitKey(10);
					  }
	
					 
		      delete  []rightdisp;
		      delete []  LRcheckresult;
		      delete []disp;
			  break;
		  }
		   case 8:
		  {
			  vector<vector<SGM_DP_TYPE>> AvgPath1_4;
			  vector<vector<SGM_DP_TYPE>> AvgPath5_8;
			  vector<vector<SGM_DP_TYPE>> AvgPath8;
			  {
			  vector<vector<SGM_DP_TYPE>> costs_L;
			  vector<vector<SGM_DP_TYPE>> costs_U;
			  vector<vector<SGM_DP_TYPE>> costs_LU;
			  vector<vector<SGM_DP_TYPE>> costs_RU;
			 
		  SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);
		  SGM_DP_U(costs_new, width, height, DMax, disp, src1, costs_U); 
		  SGM_DP_LU(costs_new, width, height, DMax, disp, src1, costs_LU);
		  SGM_DP_RU(costs_new, width, height, DMax, disp, src1, costs_RU);
		    SGM_DP_avg4(costs_L, costs_U, costs_LU,costs_RU,AvgPath1_4, width, height, DMax);
			  }

			   {
			  vector<vector<SGM_DP_TYPE>> costs_R;
			  vector<vector<SGM_DP_TYPE>> costs_D;
			  vector<vector<SGM_DP_TYPE>> costs_LD;
			  vector<vector<SGM_DP_TYPE>> costs_RD; 
		    SGM_DP_R(costs_new, width, height, DMax, disp, src1, costs_R);
		    SGM_DP_D(costs_new, width, height, DMax, disp, src1, costs_D); 
		    SGM_DP_LD(costs_new, width, height, DMax, disp, src1, costs_LD);
		    SGM_DP_RD(costs_new, width, height, DMax, disp, src1, costs_RD);
		    SGM_DP_avg4(costs_R, costs_D, costs_LD,costs_RD,AvgPath5_8, width, height, DMax);
			  }
		 SGM_DP_avg2(AvgPath1_4, AvgPath5_8, AvgPath8,width, height, DMax);
			   
		 
		 
		
		  GetDisp(AvgPath8, width, height, DMax, disp);
		   rawdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(disp, rawdisp , height, width);
		//先得到右视差图
			unsigned int *rightdisp=new unsigned int [width*height];
			 GetRightDisp(rightdisp,AvgPath8,  width,  height,DMax);
			 
			
		//再LRcheck
			unsigned int *LRcheckresult=new unsigned int [width*height];
			LRcheck(height,  width, disp, rightdisp, LRcheckresult, DMax);
			LRcheckdisp=Mat::zeros(height,width,CV_8UC1);
			 uintToMat_8UC1(LRcheckresult, LRcheckdisp , height, width);
			 //亚像素插值
			
	        subPixelrefine(LRcheckdisp,AvgPath8,subrefineout);
		//中值滤波

			medianBlur(subrefineout,finaldisp,3);

					  if(debug_display)
					  {
						cv::Mat colorresult = cv::Mat::zeros(height, width, CV_8UC3);//为彩色输出申请内存
						cv::Mat tmp=LRcheckdisp.clone();
						medianBlur(LRcheckdisp,tmp,3);
						for(int row=0;row<height;row++)
							for(int col=0;col<width;col++)
							{
							tmp.at<uchar>(row,col)*=(float)256/DMax;//把视差值线性映射到0-255，只为显示。
							}
						  GenerateFalseMap(tmp,colorresult);
						  imshow(" ",colorresult);
						  imwrite("colordisp.png",colorresult);
					//  ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
					 //  ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
					   cvWaitKey(10);
					  }
	
					
		      delete  []rightdisp;
		      delete []  LRcheckresult;
		      delete []disp;
			  break;
		  }
	  default :
		  {
		  cout<< "check the input path, path must be 1, 2 , 4 ,5 or 8. "<<endl;
		  break;

		  }
	  
	  
	  }
	
	return 0;
	}
int SGM::Getdisprity_float(float*left, float*right, int path, int DMax, bool debug_display, int width, int height, float*lrcheckdisp, float*rawdisp)
{

	uchar *src1 = new uchar[width*height];
	uchar *src2 = new uchar[width*height];
	uchar *ptl = src1;
	uchar *ptr = src2;
	for (int row = 0; row < height; row++)
		for (int col = 0; col < width; col++)
		{
			*(ptl++) = *(left++);
			*(ptr++) = *(right++);
		}

	unsigned int *ct1 = new unsigned int[width*height];
	unsigned int *ct2 = new unsigned int[width*height];
	vector< vector<int> > cpds;

	CalulateCT(src1, width, height, ct1);
	CalulateCT(src2, width, height, ct2);

	hamin(ct1, ct2, DMax, width, height, 100, cpds);
	//cout << " --> hamming cost" << endl;

	vector<vector<float>> costs_new;
	for (int i = 0; i < cpds.size(); i++)
	{
		vector<SGM_DP_TYPE> cost_new;
		for (int u = 0; u < DMax; u++)
		{

			cost_new.push_back(cpds[i][u]);
		}
		costs_new.push_back(cost_new);
	}


	unsigned int *disp = new unsigned int[width*height];

	switch (path)
	{
	case 1:
	{
		vector<vector<SGM_DP_TYPE>> costs_L;
		SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);


		GetDisp(costs_L, width, height, DMax, disp);


		uintTofloat(disp, rawdisp, height, width);
		//先得到右视差图
		unsigned int *rightdisp = new unsigned int[width*height];
		GetRightDisp(rightdisp, costs_L, width, height, DMax);


		//再LRcheck
		unsigned int *LRcheckresult = new unsigned int[width*height];
		LRcheck(height, width, disp, rightdisp, LRcheckresult, DMax);

		uintTofloat(LRcheckresult, lrcheckdisp, height, width);

		//if (debug_display)
		//{
		//	ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
		//	ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
		//	cvWaitKey(10);
		//}
		delete[]rightdisp;
		delete[]  LRcheckresult;
		delete[]disp;
		break;
	}
	case 2:
	{
		vector<vector<SGM_DP_TYPE>> costs_L;
		vector<vector<SGM_DP_TYPE>> costs_U;
		vector<vector<SGM_DP_TYPE>> costs_LU;
		SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);
		SGM_DP_U(costs_new, width, height, DMax, disp, src1, costs_U);
		SGM_DP_avg2(costs_L, costs_U, costs_LU, width, height, DMax);
		GetDisp(costs_LU, width, height, DMax, disp);
		uintTofloat(disp, rawdisp, height, width);

		//先得到右视差图
		unsigned int *rightdisp = new unsigned int[width*height];
		GetRightDisp(rightdisp, costs_LU, width, height, DMax);


		//再LRcheck
		unsigned int *LRcheckresult = new unsigned int[width*height];
		LRcheck(height, width, disp, rightdisp, LRcheckresult, DMax);
		uintTofloat(LRcheckresult, lrcheckdisp, height, width);
		Mat show = Mat::zeros(height, width, CV_8UC1);

		//if (debug_display)
		//{
		//	ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
		//	ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
		//	cvWaitKey(10);
		//}
		delete[]rightdisp;
		delete[]  LRcheckresult;
		delete[]disp;
		break;
	}
	case 4:
	{
		vector<vector<SGM_DP_TYPE>> costs_L;
		vector<vector<SGM_DP_TYPE>> costs_U;
		vector<vector<SGM_DP_TYPE>> costs_LU;
		vector<vector<SGM_DP_TYPE>> costs_RU;
		vector<vector<SGM_DP_TYPE>> AvgPath4;
		SGM_DP_L(costs_new, width, height, DMax, disp, src1, costs_L);
		SGM_DP_U(costs_new, width, height, DMax, disp, src1, costs_U);
		SGM_DP_LU(costs_new, width, height, DMax, disp, src1, costs_LU);
		SGM_DP_RU(costs_new, width, height, DMax, disp, src1, costs_RU);
		SGM_DP_avg4(costs_L, costs_U, costs_LU, costs_RU, AvgPath4, width, height, DMax);
		GetDisp(AvgPath4, width, height, DMax, disp);
		uintTofloat(disp, rawdisp, height, width);
		//先得到右视差图
		unsigned int *rightdisp = new unsigned int[width*height];
		GetRightDisp(rightdisp, AvgPath4, width, height, DMax);


		//再LRcheck
		unsigned int *LRcheckresult = new unsigned int[width*height];
		LRcheck(height, width, disp, rightdisp, LRcheckresult, DMax);

		uintTofloat(LRcheckresult, lrcheckdisp, height, width);

		//if (debug_display)
		//{
		//	ShowcolorDisp(disp, width, height, "disp", "disp.jpg", DMax);
		//	ShowcolorDisp(LRcheckresult, width, height, "LRcheckdisp", "LRcheckdisp.jpg", DMax);
		//	cvWaitKey(10);
		//}
		delete[]rightdisp;
		delete[]  LRcheckresult;
		delete[]disp;
		break;
	}
	default:
	{
		cout << "check the input path, path must be 1, 2 or 4 . " << endl;
		break;
	}


	}

	return 0;
}
#include "SGM.h"
#include "pointcloud.h"
#include "space_intersection.h"

double angleToRadian(double angle)
{
	return angle / 180.0*CV_PI;
}

int main()
{
	cv::Mat img_l = cv::imread("l.png");
	cv::Mat img_r = cv::imread("r.png");

	if (img_l.empty() || img_r.empty())
	{
		std::cout << "image path error." << std::endl;
		system("pause");

		return -1;
	}

	cv::Mat gray_l, gray_r;
	cv::cvtColor(img_l, gray_l, CV_BGR2GRAY);
	cv::cvtColor(img_r, gray_r, CV_BGR2GRAY);

	cv::resize(gray_l, gray_l, cv::Size(img_l.cols / 2, img_l.rows / 2));
	cv::resize(gray_r, gray_r, cv::Size(img_r.cols / 2, img_r.rows / 2));

	std::cout << "SGM立体匹配中.." << std::endl;
	SGM sgm;
	cv::Mat disp;
	sgm.GetDisprity_mat(gray_l, gray_r, 4, 256, false, disp);
	cv::resize(disp, disp, cv::Size(disp.cols * 2, disp.rows * 2));

	cv::medianBlur(disp, disp, 7);

	std::vector<SpaceIntersection::Matcher>matches;
	for (auto i = 0; i < disp.rows; ++i)
	{
		float *data_src = disp.ptr<float>(i);
		for (auto j = 0; j < disp.cols; ++j)
		{
			if (data_src[j] > 0)
			{
				data_src[j] = 2 * data_src[j];

				SpaceIntersection::Matcher m(cv::Point2f(j, i), cv::Point2f(fabs(j - data_src[j]), i));
				matches.push_back(m);
			}
		}
	}
	std::cout << "SGM立体匹配已完成！" << std::endl;

	/*
	*
	* 前方交会部分
	*
	*/
	std::cout << "开始前方交会计算点云.." << std::endl;
	SpaceIntersection insec(matches);
	insec.setInteriorOrientation(4166.67, 4166.67, 4166.67, 0, 0);
	insec.setExteriorOrientation_l(2240.000000, 2002.000000, 8000.000000, 0, 0, 0).setExteriorOrientation_r(2740.000000, 2002.000000, 8000.000000, 0, 0, 0);
	insec.setImgSize(img_l.cols, img_l.rows);

	insec.compute(img_l).storage("点云.txt");
	std::cout << "前方交会计算点云已完成！" << std::endl;

	system("pause");
	return 0;
}
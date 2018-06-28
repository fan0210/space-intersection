#ifndef SPACE_INTERSECTION_H_
#define SPACE_INTERSECTION_H_

#include <opencv.hpp>

#include "pointcloud.h"

//单个影像对前方交会生成点云
class SpaceIntersection
{
public:
	class Matcher
	{
	public:
		using Point = cv::Point2f;
		Matcher() = default;
		Matcher(const Point &pt_l, const Point &pt_r) :m_ptl(pt_l), m_ptr(pt_r) {}

		const Point &getPt_l()const { return m_ptl; }
		const Point &getPt_r()const { return m_ptr; }
	private:
		Point m_ptl;   //左影像上的点
		Point m_ptr;   //右影像上的点
	};

	using Point = PointCloud::Point;

	SpaceIntersection(std::vector<Matcher> &matches) :m_matches(matches) {}

	SpaceIntersection &setInteriorOrientation(double f, double fx, double fy, double x0, double y0);                  //设置内方位元素
	SpaceIntersection &setExteriorOrientation_l(double X_l, double Y_l, double Z_l, double phi_l, double omega_l, double kappa_l);        //设置左影像外方位元素
	SpaceIntersection &setExteriorOrientation_r(double X_r, double Y_r, double Z_r, double phi_r, double omega_r, double kappa_r);        //设置右影像外方位元素
	SpaceIntersection &setImgSize(int imgWidth, int imgHeight);           //设置图像的大小         

	PointCloud &compute(const cv::Mat &rgbImg = cv::Mat());
	const PointCloud &getPointCloud()const { return m_pd; }
private:
	std::vector<Matcher> &m_matches;
	PointCloud m_pd;

	double m_f = 0, m_fx = 0, m_fy = 0, m_x0 = 0, m_y0 = 0;
	double m_Xl, m_Yl, m_Zl, m_phil, m_omegal, m_kappal;
	double m_Xr, m_Yr, m_Zr, m_phir, m_omegar, m_kappar;
	int m_width = 0, m_height = 0;

	double R_l[9]{ 0 };
	double R_r[9]{ 0 };

	Point singleMatchCompute(const Matcher &, const cv::Mat &)const;
	
	void correctXY(const cv::Point2f &, cv::Point2f &)const;         //改正像素坐标系下的坐标到以像素为单位的像平面坐标系下
	
	void transform_l(const Point &, Point &)const;             //左像空间坐标系下的点转换到像空间辅助坐标系下
	void transform_r(const Point &, Point &)const;             //右像空间坐标系下的点转换到像空间辅助坐标系下

	Point projectionCoefficientSlove(const Matcher &, const cv::Mat &)const;              //点投影系数法
	Point leastSquaresSlove(const Matcher &, const cv::Mat &)const;                  //最小二乘迭代求解精确三维坐标
};  
#endif

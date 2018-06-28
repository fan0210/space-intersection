#include "space_intersection.h"

void SpaceIntersection::correctXY(const cv::Point2f &inPoint, cv::Point2f &outPoint)const
{
	double x0 = m_x0 + (m_width - 1) / 2.0;
	double y0 = m_y0 + (m_height - 1) / 2.0;

	outPoint.x = (inPoint.x - x0)*m_fx / m_f;
	outPoint.y = (inPoint.y - y0)*m_fy / m_f;
}

void SpaceIntersection::transform_l(const Point &inPoint, Point &outPoint)const
{
	double x = R_l[0] * inPoint.getX() + R_l[1] * inPoint.getY() + R_l[2] * inPoint.getZ();
	double y = R_l[3] * inPoint.getX() + R_l[4] * inPoint.getY() + R_l[5] * inPoint.getZ();
	double z = R_l[6] * inPoint.getX() + R_l[7] * inPoint.getY() + R_l[8] * inPoint.getZ();

	outPoint = Point(x, y, z);
}

void SpaceIntersection::transform_r(const Point &inPoint, Point &outPoint)const
{
	double x = R_r[0] * inPoint.getX() + R_r[1] * inPoint.getY() + R_r[2] * inPoint.getZ();
	double y = R_r[3] * inPoint.getX() + R_r[4] * inPoint.getY() + R_r[5] * inPoint.getZ();
	double z = R_r[6] * inPoint.getX() + R_r[7] * inPoint.getY() + R_r[8] * inPoint.getZ();

	outPoint = Point(x, y, z);
}

SpaceIntersection::Point SpaceIntersection::projectionCoefficientSlove(const Matcher &m, const cv::Mat &rgbImg)const
{
	bool rgb = !rgbImg.empty();

	cv::Point2f outPoint_l, outPoint_r;
	cv::Point2f inPoint_l(m.getPt_l().x, m.getPt_l().y);
	cv::Point2f inPoint_r(m.getPt_r().x, m.getPt_r().y);

	correctXY(inPoint_l, outPoint_l);
	correctXY(inPoint_r, outPoint_r);

	Point pt_l(outPoint_l.x, outPoint_l.y, -m_f);
	Point pt_r(outPoint_r.x, outPoint_r.y, -m_f);

	Point pt_l_, pt_r_;
	transform_l(pt_l, pt_l_);
	transform_r(pt_r, pt_r_);

	double Bx = m_Xr - m_Xl;
	double By = m_Yr - m_Yl;
	double Bz = m_Zr - m_Zl;
	
	double N1_Y = (Bx*pt_r_.getZ() - Bz*pt_r_.getX()) / (pt_l_.getX()*pt_r_.getZ() - pt_r_.getX()*pt_l_.getZ());
	double N2_Y = (Bx*pt_l_.getZ() - Bz*pt_l_.getX()) / (pt_l_.getX()*pt_r_.getZ() - pt_r_.getX()*pt_l_.getZ());
	double y1_Y = m_Yl + N1_Y*pt_l_.getY();
	double y2_Y = m_Yr + N2_Y*pt_r_.getY();
	double x_Y = m_Xl + N1_Y*pt_l_.getX();
	double z_Y = m_Zl + N1_Y*pt_l_.getZ();

	double 	x = x_Y;
	double 	z = z_Y;
	double 	y = (y1_Y + y2_Y) / 2.0;

	if (!rgb)
	{
		unsigned char r = rgbImg.at<cv::Vec3b>(int(m.getPt_l().y), int(m.getPt_l().x))[2];
		unsigned char g = rgbImg.at<cv::Vec3b>(int(m.getPt_l().y), int(m.getPt_l().x))[1];
		unsigned char b = rgbImg.at<cv::Vec3b>(int(m.getPt_l().y), int(m.getPt_l().x))[0];

		return Point(x, y, z, r, g, b);
	}
	else
		return Point(x, y, z);
}

SpaceIntersection::Point SpaceIntersection::leastSquaresSlove(const Matcher &m, const cv::Mat &rgbImg)const
{
	bool rgb = !rgbImg.empty();

	double X[3]{ 0 };
	cv::Mat_<double> A(cv::Size(3, 4), 0);
	cv::Mat_<double> L(cv::Size(1, 4), 0);
	size_t iteratorTimes = 0;

	Point pt_rough = projectionCoefficientSlove(m, rgbImg);          //由点投影系数法获取初值
	if (pt_rough.getX() != 0 || pt_rough.getX() != 0 || pt_rough.getX() != 0)
	{
		X[0] = pt_rough.getX();
		X[1] = pt_rough.getY();
		X[2] = pt_rough.getZ();

		cv::Point2f outPoint_l, outPoint_r;
		cv::Point2f inPoint_l(m.getPt_l().x, m.getPt_l().y);
		cv::Point2f inPoint_r(m.getPt_r().x, m.getPt_r().y);

		correctXY(inPoint_l, outPoint_l);
		correctXY(inPoint_r, outPoint_r);

		cv::Mat_<double> dX_(cv::Size(1, 3), 100000000);
		while (iteratorTimes < 5)
		{
			double X_l = (R_l[0] * (X[0] - m_Xl) + R_l[3] * (X[1] - m_Yl) + R_l[6] * (X[2] - m_Zl));
			double Y_l = (R_l[1] * (X[0] - m_Xl) + R_l[4] * (X[1] - m_Yl) + R_l[7] * (X[2] - m_Zl));
			double Z_l = (R_l[2] * (X[0] - m_Xl) + R_l[5] * (X[1] - m_Yl) + R_l[8] * (X[2] - m_Zl));
			double X_r = (R_r[0] * (X[0] - m_Xr) + R_r[3] * (X[1] - m_Yr) + R_r[6] * (X[2] - m_Zr));
			double Y_r = (R_r[1] * (X[0] - m_Xr) + R_r[4] * (X[1] - m_Yr) + R_r[7] * (X[2] - m_Zr));
			double Z_r = (R_r[2] * (X[0] - m_Xr) + R_r[5] * (X[1] - m_Yr) + R_r[8] * (X[2] - m_Zr));

			A.at<double>(0, 0) = -(R_l[0] * m_f + R_l[2] * outPoint_l.x) / Z_l;
			A.at<double>(0, 1) = -(R_l[3] * m_f + R_l[5] * outPoint_l.x) / Z_l;
			A.at<double>(0, 2) = -(R_l[6] * m_f + R_l[8] * outPoint_l.x) / Z_l;
			A.at<double>(1, 0) = -(R_l[1] * m_f + R_l[2] * outPoint_l.y) / Z_l;
			A.at<double>(1, 1) = -(R_l[4] * m_f + R_l[5] * outPoint_l.y) / Z_l;
			A.at<double>(1, 2) = -(R_l[7] * m_f + R_l[8] * outPoint_l.y) / Z_l;
			A.at<double>(2, 0) = -(R_r[0] * m_f + R_r[2] * outPoint_r.x) / Z_r;
			A.at<double>(2, 1) = -(R_r[3] * m_f + R_r[5] * outPoint_r.x) / Z_r;
			A.at<double>(2, 2) = -(R_r[6] * m_f + R_r[8] * outPoint_r.x) / Z_r;
			A.at<double>(3, 0) = -(R_r[1] * m_f + R_r[2] * outPoint_r.y) / Z_r;
			A.at<double>(3, 1) = -(R_r[4] * m_f + R_r[5] * outPoint_r.y) / Z_r;
			A.at<double>(3, 2) = -(R_r[7] * m_f + R_r[8] * outPoint_r.y) / Z_r;

			double x0_l = -m_f*X_l / Z_l;
			double y0_l = -m_f*Y_l / Z_l;
			double x0_r = -m_f*X_r / Z_r;
			double y0_r = -m_f*Y_r / Z_r;

			L.at<double>(0, 0) = outPoint_l.x - x0_l;
			L.at<double>(1, 0) = outPoint_l.y - y0_l;
			L.at<double>(2, 0) = outPoint_r.x - x0_r;
			L.at<double>(3, 0) = outPoint_r.y - y0_r;

			cv::Mat_<double> dX = (A.t()*A).inv()*(A.t()*L);
			if (fabs(dX.at<double>(0, 0) - dX_.at<double>(0, 0)) < 0.01
				&&fabs(dX.at<double>(1, 0) - dX_.at<double>(1, 0)) < 0.01
				&&fabs(dX.at<double>(2, 0) - dX_.at<double>(2, 0)) < 0.01)
			{
				X[0] += dX.at<double>(0, 0);
				X[1] += dX.at<double>(1, 0);
				X[2] += dX.at<double>(2, 0);

				if (rgb)
				{
					unsigned char r = rgbImg.at<cv::Vec3b>(int(m.getPt_l().y), int(m.getPt_l().x))[2];
					unsigned char g = rgbImg.at<cv::Vec3b>(int(m.getPt_l().y), int(m.getPt_l().x))[1];
					unsigned char b = rgbImg.at<cv::Vec3b>(int(m.getPt_l().y), int(m.getPt_l().x))[0];

					return Point(X[0], X[1], X[2], r, g, b);
				}
				else
					return Point(X[0], X[1], X[2]);
			}
			else
			{
				X[0] += dX.at<double>(0, 0);
				X[1] += dX.at<double>(1, 0);
				X[2] += dX.at<double>(2, 0);

				dX_ = dX;
			}

			++iteratorTimes;
		}
	}
	return Point(0, 0, 0);
}

SpaceIntersection &SpaceIntersection::setInteriorOrientation(double f, double fx, double fy, double x0, double y0)
{
	m_f = f;
	m_fx = fx;
	m_fy = fy;
	m_x0 = x0;
	m_y0 = y0;

	return *this;
}

SpaceIntersection &SpaceIntersection::setExteriorOrientation_l(double X_l, double Y_l, double Z_l, double phi_l, double omega_l, double kappa_l)
{
	m_Xl = X_l;
	m_Yl = Y_l;
	m_Zl = Z_l;
	m_phil = phi_l;
	m_omegal = omega_l;
	m_kappal = kappa_l;

	return *this;
}

SpaceIntersection &SpaceIntersection::setExteriorOrientation_r(double X_r, double Y_r, double Z_r, double phi_r, double omega_r, double kappa_r)
{
	m_Xr = X_r;
	m_Yr = Y_r;
	m_Zr = Z_r;
	m_phir = phi_r;
	m_omegar = omega_r;
	m_kappar = kappa_r;

	return *this;
}

SpaceIntersection &SpaceIntersection::setImgSize(int imgWidth, int imgHeight)
{
	m_width = imgWidth;
	m_height = imgHeight;

	return *this;
}

inline SpaceIntersection::Point SpaceIntersection::singleMatchCompute(const Matcher &m,const cv::Mat &rgbImg)const
{
	return leastSquaresSlove(m, rgbImg);
}

PointCloud &SpaceIntersection::compute(const cv::Mat &rgbImg)
{
	R_l[0] = cos(m_phil)*cos(m_kappal) - sin(m_phil)*sin(m_omegal)*sin(m_kappal);
	R_l[1] = -cos(m_phil)*sin(m_kappal) - sin(m_phil)*sin(m_omegal)*cos(m_kappal);
	R_l[2] = -sin(m_phil)*cos(m_omegal);
	R_l[3] = cos(m_omegal)*sin(m_kappal);
	R_l[4] = cos(m_omegal)*cos(m_kappal);
	R_l[5] = -sin(m_omegal);
	R_l[6] = sin(m_phil)*cos(m_kappal) + cos(m_phil)*sin(m_omegal)*sin(m_kappal);
	R_l[7] = -sin(m_phil)*sin(m_kappal) + cos(m_phil)*sin(m_omegal)*cos(m_kappal);
	R_l[8] = cos(m_phil)*cos(m_omegal);

	R_r[0] = cos(m_phir)*cos(m_kappar) - sin(m_phir)*sin(m_omegar)*sin(m_kappar);
	R_r[1] = -cos(m_phir)*sin(m_kappar) - sin(m_phir)*sin(m_omegar)*cos(m_kappar);
	R_r[2] = -sin(m_phir)*cos(m_omegar);
	R_r[3] = cos(m_omegar)*sin(m_kappar);
	R_r[4] = cos(m_omegar)*cos(m_kappar);
	R_r[5] = -sin(m_omegar);
	R_r[6] = sin(m_phir)*cos(m_kappar) + cos(m_phir)*sin(m_omegar)*sin(m_kappar);
	R_r[7] = -sin(m_phir)*sin(m_kappar) + cos(m_phir)*sin(m_omegar)*cos(m_kappar);
	R_r[8] = cos(m_phir)*cos(m_omegar);

	for (auto it = m_matches.cbegin(); it != m_matches.cend(); ++it)
	{
		Point pt = singleMatchCompute(*it, rgbImg);
		if (pt.getX() != 0 || pt.getY() != 0 || pt.getZ() != 0)
			m_pd.push_back(pt);
	}

	return m_pd;
}
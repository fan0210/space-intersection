#ifndef POINT_CLOUD_H_
#define POINT_CLOUD_H_

#include <vector>
#include <memory>
#include <fstream>

class PointCloud
{
public:
	class Point
	{
	public:
		Point() = default;
		Point(float x, float y, float z, unsigned char r = 255, unsigned char g = 255, unsigned char b = 255) :m_x(x), m_y(y), m_z(z), m_r(r), m_g(g), m_b(b) {}

		float getX()const { return m_x; }
		float getY()const { return m_y; }
		float getZ()const { return m_z; }

		int getR()const { return m_r; }
		int getG()const { return m_g; }
		int getB()const { return m_b; }
	private:
		float m_x = 0;
		float m_y = 0;
		float m_z = 0;

		unsigned char m_r = 0;
		unsigned char m_g = 0;
		unsigned char m_b = 0;
	};

	PointCloud() = default;
	explicit PointCloud(const std::vector<Point> &pts) :m_ptsPtr(std::make_shared<std::vector<Point>>(pts)) {}

	size_t size()const { return m_ptsPtr->size(); }
	bool empty()const { return m_ptsPtr->empty(); }
	void clear() { m_ptsPtr->clear(); }

	auto &begin() { return m_ptsPtr->begin(); }
	auto &end() { return m_ptsPtr->end(); }

	auto &cbegin()const { return m_ptsPtr->cbegin(); }
	auto &cend()const { return m_ptsPtr->cend(); }

	const std::vector<Point> &getPts()const;
	void push_back(const Point &pt);
	const PointCloud clone()const;
	void storage(const std::string &fileAddr)const;
	void storage(std::fstream &fout)const;

	std::shared_ptr<std::vector<Point>> &data() { return m_ptsPtr; }
	const std::shared_ptr<std::vector<Point>> &cdata()const { return m_ptsPtr; }
private:
	std::shared_ptr<std::vector<Point>> m_ptsPtr = std::make_shared<std::vector<Point>>();
};

#endif

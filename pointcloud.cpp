#include "pointcloud.h"

const std::vector<PointCloud::Point> &PointCloud::getPts()const
{
	return *m_ptsPtr;
}

void PointCloud::push_back(const Point &pt)
{
	m_ptsPtr->push_back(pt);
}

const PointCloud PointCloud::clone()const
{
	PointCloud pt(*m_ptsPtr);

	return pt;
}

void PointCloud::storage(const std::string &fileAddr)const
{
	std::fstream file(fileAddr, std::ios::app);
	if (!file)
		return;

	for (auto it = m_ptsPtr->begin(); it != m_ptsPtr->end(); ++it)
		file << it->getX() << "   " << it->getY() << "   " << it->getZ() << "   " << it->getR() << "   " << it->getG() << "   " << it->getB() << std::endl;
}

void PointCloud::storage(std::fstream &fout)const
{
	for (auto it = m_ptsPtr->begin(); it != m_ptsPtr->end(); ++it)
		fout << it->getX() << "   " << it->getY() << "   " << it->getZ() << "   " << it->getR() << "   " << it->getG() << "   " << it->getB() << std::endl;
}
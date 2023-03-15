/*
Copyright 2011. All rights reserved.
Institute of Measurement and Control Systems
Karlsruhe Institute of Technology, Germany

Authors: Andreas Geiger

libicp is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or any later version.

libicp is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
libicp; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301, USA
*/

#include "icp.h"
using namespace std;

Icp::Icp(const std::vector<Eigen::Vector3f>& M)
    :m_max_iter(200), m_min_delta(1e-4)
{

    // check for minimum number of points
    if (M.size() < 5) {
        cout << "ERROR: LIBICP works only with at least 5 model points" << endl;
        m_kd_tree = 0;
        return;
    }

    // copy model points to M_data
    m_kd_data.resize(boost::extents[M.size()][2]);
    for (int32_t m = 0; m < M.size(); m++)
    {
        //for (int32_t n = 0; n < 2; n++)
        m_kd_data[m][0] = M[m].x();
        m_kd_data[m][1] = M[m].y();
    }

    // build a kd tree from the model point cloud
    m_kd_tree = new kdtree::KDTree(m_kd_data);
}

Icp::~Icp() {
    if (m_kd_tree)
        delete m_kd_tree;
}

double Icp::fit(const std::vector<Eigen::Vector2f>& T, Eigen::Isometry2f& R, float indist = -1.0f)
{
    // make sure we have a model tree
    if (!m_kd_tree) {
        cout << "ERROR: No model available." << endl;
        return 0;
    }

    // check for minimum number of points
    if (T_num < 5) {
        cout << "ERROR: Icp works only with at least 5 template points" << endl;
        return 0;
    }

    // set active points
    vector<int32_t> active;
    if (indist <= 0) {
        active.clear();
        for (int32_t i = 0; i < T_num; i++)
            active.push_back(i);
    }
    else {
        active = getInliers(T, T_num, R, t, indist);
    }

    // run icp
    fitIterate(T, T_num, R, t, indist);

    return getResidual(T, T_num, R, t, m_active);
}

void Icp::fitIterate(double* T, const int32_t T_num, Matrix& R, Matrix& t, double indist /*= -1*/)
{
    if (indist <= 0) {
        m_active.clear(); m_active.resize(T_num);
        for (int32_t i = 0; i < T_num; i++) {
            m_active[i] = i;
        }
        m_inlier_ratio = 1;
    }
    double delta = 1000;
    int32_t iter;
    for (iter = 0; iter<m_max_iter && delta>m_min_delta; iter++) {
        if (indist > 0) {
            indist = std::max(indist * 0.9, 0.05);
            m_active = getInliers(T, T_num, R, t, indist);
            m_inlier_ratio = (double)m_active.size() / T_num;
        }
        delta = fitStep(T, T_num, R, t, m_active);
    }
}

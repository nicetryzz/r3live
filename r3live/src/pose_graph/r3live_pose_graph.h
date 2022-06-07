#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Core>


#include "rgb_map/image_frame.hpp"
#include "rgb_map/pointcloud_rgbd.hpp"
#include "rgb_map/offline_map_recorder.hpp"

#include <g2o/core/base_vertex.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>

#include <sophus/se3.hpp>

typedef Eigen::Matrix<double, 6, 6> Matrix6d;

struct PoseGraphVertex{
    Eigen::Matrix3d R;
    Eigen::Vector3d C;
} typedef PoseGraphVertex;

struct PoseGraphEdge{
    int id1,id2;
    Eigen::Matrix3d R;
    Eigen::Vector3d C;
} typedef PoseGraphEdge;

bool pose_graph( Offline_map_recorder &r3live_map_recorder);
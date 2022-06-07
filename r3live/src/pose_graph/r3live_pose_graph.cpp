#include "r3live_pose_graph.h"

using namespace std;
using namespace Eigen;
using Sophus::SE3d;
using Sophus::SO3d;

/************************************************
 * 本程序演示如何用g2o solver进行位姿图优化
 * sphere.g2o是人工生成的一个Pose graph，我们来优化它。
 * 尽管可以直接通过load函数读取整个图，但我们还是自己来实现读取代码，以期获得更深刻的理解
 * 本节使用李代数表达位姿图，节点和边的方式为自定义
 * **********************************************/

// 给定误差求J_R^{-1}的近似
Matrix6d JRInv(const SE3d &e) {
    Matrix6d J;
    J.block(0, 0, 3, 3) = SO3d::hat(e.so3().log());
    J.block(0, 3, 3, 3) = SO3d::hat(e.translation());
    J.block(3, 0, 3, 3) = Matrix3d::Zero(3, 3);
    J.block(3, 3, 3, 3) = SO3d::hat(e.so3().log());
    // J = J * 0.5 + Matrix6d::Identity();
    J = Matrix6d::Identity();    // try Identity if you want
    return J;
}

// 李代数顶点
typedef Matrix<double, 6, 1> Vector6d;

class VertexSE3LieAlgebra : public g2o::BaseVertex<6, SE3d> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(istream &is) override {
        double data[7];
        for (int i = 0; i < 7; i++)
            is >> data[i];
        setEstimate(SE3d(
            Quaterniond(data[6], data[3], data[4], data[5]),
            Vector3d(data[0], data[1], data[2])
        ));
        return true;
    }

    virtual bool write(ostream &os) const override {
        os << id() << " ";
        Quaterniond q = _estimate.unit_quaternion();
        os << _estimate.translation().transpose() << " ";
        os << q.coeffs()[0] << " " << q.coeffs()[1] << " " << q.coeffs()[2] << " " << q.coeffs()[3] << endl;
        return true;
    }
    
    bool setVertex(PoseGraphVertex v)  {
        Quaterniond q = Quaterniond(v.R);
        q.normalize();
        setEstimate(SE3d(q, v.C));
        return true;
    }

    PoseGraphVertex getVertex() {
        PoseGraphVertex v;
        Quaterniond q = _estimate.unit_quaternion();
        v.R = q.toRotationMatrix();
        v.C = _estimate.translation().transpose();
        return v;
    }

    virtual void setToOriginImpl() override {
        _estimate = SE3d();
    }

    // 左乘更新
    virtual void oplusImpl(const double *update) override {
        Vector6d upd;
        upd << update[0], update[1], update[2], update[3], update[4], update[5];
        _estimate = SE3d::exp(upd) * _estimate;
    }
};

// 两个李代数节点之边
class EdgeSE3LieAlgebra : public g2o::BaseBinaryEdge<6, SE3d, VertexSE3LieAlgebra, VertexSE3LieAlgebra> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    virtual bool read(istream &is) override {
        double data[7];
        for (int i = 0; i < 7; i++)
            is >> data[i];
        Quaterniond q(data[6], data[3], data[4], data[5]);
        q.normalize();
        setMeasurement(SE3d(q, Vector3d(data[0], data[1], data[2])));
        for (int i = 0; i < information().rows() && is.good(); i++)
            for (int j = i; j < information().cols() && is.good(); j++) {
                is >> information()(i, j);
                if (i != j)
                    information()(j, i) = information()(i, j);
            }
        return true;
    }

    virtual bool write(ostream &os) const override {
        VertexSE3LieAlgebra *v1 = static_cast<VertexSE3LieAlgebra *> (_vertices[0]);
        VertexSE3LieAlgebra *v2 = static_cast<VertexSE3LieAlgebra *> (_vertices[1]);
        os << v1->id() << " " << v2->id() << " ";
        SE3d m = _measurement;
        Quaterniond q = m.unit_quaternion();
        os << m.translation().transpose() << " ";
        os << q.coeffs()[0] << " " << q.coeffs()[1] << " " << q.coeffs()[2] << " " << q.coeffs()[3] << " ";

        // information matrix 
        for (int i = 0; i < information().rows(); i++)
            for (int j = i; j < information().cols(); j++) {
                os << information()(i, j) << " ";
            }
        os << endl;
        return true;
    }

    bool setEdge(PoseGraphEdge e,int value) {
        Quaterniond q = Quaterniond(e.R);
        q.normalize();
        setMeasurement(SE3d(q, e.C));
        for (int i = 0; i < information().rows(); i++){
            for (int j = i; j < information().cols(); j++)  {
                if(i == j && i < 3)
                information()(i,j) = value;
                else if(i == j)
                information()(i,j) = value;
                else{
                    information()(j, i) = 0;
                    information()(i, j) = 0;
                }
            }
        }
        return true;
    }

    PoseGraphEdge getEdge() {
        PoseGraphEdge e;
        VertexSE3LieAlgebra *v1 = static_cast<VertexSE3LieAlgebra *> (_vertices[0]);
        VertexSE3LieAlgebra *v2 = static_cast<VertexSE3LieAlgebra *> (_vertices[1]);
        e.id1 = v1->id();
        e.id2 = v2->id();
        SE3d m = _measurement;
        Quaterniond q = m.unit_quaternion();
        e.R = q.toRotationMatrix();
        e.C = m.translation().transpose();
        return e;
    }

    // 误差计算与书中推导一致
    virtual void computeError() override {
        SE3d v1 = (static_cast<VertexSE3LieAlgebra *> (_vertices[0]))->estimate();
        SE3d v2 = (static_cast<VertexSE3LieAlgebra *> (_vertices[1]))->estimate();
        _error = (_measurement.inverse() * v1.inverse() * v2).log();
    }

    // 雅可比计算
    virtual void linearizeOplus() override {
        SE3d v1 = (static_cast<VertexSE3LieAlgebra *> (_vertices[0]))->estimate();
        SE3d v2 = (static_cast<VertexSE3LieAlgebra *> (_vertices[1]))->estimate();
        Matrix6d J = JRInv(SE3d::exp(_error));
        // 尝试把J近似为I？
        _jacobianOplusXi = -J * v2.inverse().Adj();
        _jacobianOplusXj = J * v2.inverse().Adj();
    }
};

bool pose_graph( Offline_map_recorder &r3live_map_recorder){

    typedef g2o::BlockSolver<g2o::BlockSolverTraits<6, 6>> BlockSolverType;
    typedef g2o::LinearSolverEigen<BlockSolverType::PoseMatrixType> LinearSolverType;
    auto solver = new g2o::OptimizationAlgorithmLevenberg(
        g2o::make_unique<BlockSolverType>(g2o::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer;     // 图模型
    optimizer.setAlgorithm(solver);   // 设置求解器
    optimizer.setVerbose(true);       // 打开调试输出

    int vertexCnt = 0, edgeCnt = 0; // 顶点和边的数量

    vector<VertexSE3LieAlgebra *> vectices;
    vector<EdgeSE3LieAlgebra *> edges;

    int     number_of_image_frame = r3live_map_recorder.m_pts_in_views_vec.size();

    for (int frame_idx = 0; frame_idx < number_of_image_frame; frame_idx++ )
    {
        std::shared_ptr< Image_frame > img_ptr = r3live_map_recorder.m_image_pose_vec[ frame_idx ];
        vec_3   pose_t = -img_ptr->m_pose_c2w_q.toRotationMatrix().transpose() * img_ptr->m_pose_c2w_t;
        VertexSE3LieAlgebra *v = new VertexSE3LieAlgebra();
        PoseGraphVertex vertex;
        vertex.R = img_ptr->m_pose_c2w_q.toRotationMatrix().transpose();
        vertex.C = pose_t;
        v->setId(frame_idx);
        v->setVertex(vertex);
        optimizer.addVertex(v);
        vectices.push_back(v);
        if(frame_idx > 0){
            EdgeSE3LieAlgebra *e = new EdgeSE3LieAlgebra();
            e->setId(frame_idx - 1);
            e->setVertex(0, optimizer.vertices()[frame_idx - 1]);
            e->setVertex(1, optimizer.vertices()[frame_idx]);
            PoseGraphEdge edge;
            VertexSE3LieAlgebra *v1 = vectices[frame_idx - 1];
            PoseGraphVertex vertex1 = v1->getVertex();
            VertexSE3LieAlgebra *v2 = vectices[frame_idx];
            PoseGraphVertex vertex2 = v2->getVertex();

            edge.R = vertex1.R.transpose()*vertex2.R;
            edge.C = vertex1.R.transpose()*(vertex2.C-vertex1.C);
            e->setEdge(edge,frame_idx);
            optimizer.addEdge(e);
            edges.push_back(e);
        }
    }
    EdgeSE3LieAlgebra *e = new EdgeSE3LieAlgebra();
    e->setId(number_of_image_frame - 1);
    e->setVertex(0, optimizer.vertices()[number_of_image_frame - 1]);
    e->setVertex(1, optimizer.vertices()[0]);
    PoseGraphEdge edge;
    edge.R = Eigen::Matrix3d::Identity();
    edge.C =  Eigen::Vector3d::Zero();
    e->setEdge(edge,number_of_image_frame);
    optimizer.addEdge(e);
    edges.push_back(e);
    optimizer.initializeOptimization();
    optimizer.optimize(30);
    PoseGraphVertex vertex_change = vectices[0]->getVertex();

    for (int frame_idx = 0; frame_idx < number_of_image_frame; frame_idx++ )
    {
        std::shared_ptr< Image_frame > img_ptr = r3live_map_recorder.m_image_pose_vec[ frame_idx ];
        VertexSE3LieAlgebra *v = vectices[frame_idx];
        PoseGraphVertex vertex = v->getVertex();
        img_ptr->m_pose_c2w_q = Quaterniond(vertex_change.R*vertex.R.transpose());
        img_ptr->m_pose_c2w_t = -(vertex_change.R*vertex.R.transpose())*(vertex.C - vertex_change.C);
    }
    return true;
}

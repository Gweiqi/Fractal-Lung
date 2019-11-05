#include <eigen/Core>

int main()
{
    Eigen::Matrix3f m = Eigen::Matrix3f::Zero();
    Eigen::Matrix3f n = Eigen::Matrix3f::Zero();
    m.block<2, 2>(0, 0) = n.block<2, 2>(0, 0);
}
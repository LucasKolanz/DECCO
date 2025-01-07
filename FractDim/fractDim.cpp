#include <open3d/geometry/PointCloud.h>
#include <open3d/geometry/KDTreeFlann.h>
#include <open3d/utility/Console.h>
#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <numeric>
#include <fstream>

using namespace open3d;
using namespace std;

// Function to sample points on a sphere surface
std::vector<Eigen::Vector3d> SamplePointsOnSphere(
    const Eigen::Vector3d& center,
    double radius,
    int num_points) {
    std::vector<Eigen::Vector3d> points;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < num_points; ++i) {
        double u = dis(gen);
        double v = dis(gen);
        double theta = 2 * M_PI * u;
        double phi = acos(2 * v - 1);
        double x = center[0] + radius * sin(phi) * cos(theta);
        double y = center[1] + radius * sin(phi) * sin(theta);
        double z = center[2] + radius * cos(phi);
        points.emplace_back(x, y, z);
    }
    return points;
}

int main() {
    // Sphere centers and radii (replace with your data)
    std::vector<Eigen::Vector3d> centers = {
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 1.0}
        // Add more centers
    };
    std::vector<double> radii = {
        0.5,
        0.3
        // Corresponding radii
    };
    int num_points_per_sphere = 1000; // Adjust as needed

    // Sample points on spheres
    std::vector<Eigen::Vector3d> all_points;

    for (size_t i = 0; i < centers.size(); ++i) {
        auto sphere_points = SamplePointsOnSphere(
            centers[i], radii[i], num_points_per_sphere);
        all_points.insert(all_points.end(), sphere_points.begin(), sphere_points.end());
    }

    // Create a PointCloud
    auto pcd = std::make_shared<geometry::PointCloud>();
    pcd->points_ = all_points;

    // Define voxel sizes
    std::vector<double> voxel_sizes;
    int num_scales = 20;
    double voxel_size_min = 0.01; // Adjust as needed
    double voxel_size_max = 1.0;  // Adjust as needed

    for (int i = 0; i < num_scales; ++i) {
        double voxel_size = voxel_size_min * pow((voxel_size_max / voxel_size_min),
                            static_cast<double>(i) / (num_scales - 1));
        voxel_sizes.push_back(voxel_size);
    }

    // Voxel counting
    std::vector<std::pair<double, int>> voxel_counts;

    for (const auto& voxel_size : voxel_sizes) {
        auto voxel_grid = geometry::VoxelGrid::CreateFromPointCloud(*pcd, voxel_size);
        int num_voxels = static_cast<int>(voxel_grid->voxels_.size());
        voxel_counts.emplace_back(voxel_size, num_voxels);
    }

    // Prepare data for linear regression
    std::vector<double> log_voxel_sizes;
    std::vector<double> log_voxel_counts;

    for (const auto& vc : voxel_counts) {
        log_voxel_sizes.push_back(log(1.0 / vc.first));
        log_voxel_counts.push_back(log(static_cast<double>(vc.second)));
    }

    // Perform linear regression
    size_t n = log_voxel_sizes.size();
    double sum_x = std::accumulate(log_voxel_sizes.begin(), log_voxel_sizes.end(), 0.0);
    double sum_y = std::accumulate(log_voxel_counts.begin(), log_voxel_counts.end(), 0.0);
    double sum_xy = std::inner_product(log_voxel_sizes.begin(), log_voxel_sizes.end(), log_voxel_counts.begin(), 0.0);
    double sum_x2 = std::inner_product(log_voxel_sizes.begin(), log_voxel_sizes.end(), log_voxel_sizes.begin(), 0.0);

    double denominator = n * sum_x2 - sum_x * sum_x;
    double slope = (n * sum_xy - sum_x * sum_y) / denominator;
    double intercept = (sum_y * sum_x2 - sum_x * sum_xy) / denominator;

    std::cout << "Estimated fractal dimension: " << slope << std::endl;

    // Write data to CSV
    std::ofstream outfile("box_counting_data.csv");
    outfile << "log(1/voxel_size),log(num_voxels)\n";
    for (size_t i = 0; i < n; ++i) {
        outfile << log_voxel_sizes[i] << "," << log_voxel_counts[i] << "\n";
    }
    outfile.close();

    return 0;
}

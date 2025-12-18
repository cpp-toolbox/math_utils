#include "math_utils.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace math_utils {

int non_neg_mod(int value, int mod) { return (value % mod + mod) % mod; }

std::pair<float, float> extract_yaw_pitch(const glm::vec3 &forward) {
    glm::vec3 dir = glm::normalize(forward);
    float pitch = std::asin(dir.y);
    float yaw = std::atan2(dir.x, dir.z);
    yaw = glm::degrees(yaw);
    pitch = glm::degrees(pitch);
    return {yaw, pitch};
}

double map_range(double value, double in_min, double in_max, double out_min, double out_max) {
    value = std::clamp(value, in_min, in_max);
    return out_min + (value - in_min) * (out_max - out_min) / (in_max - in_min);
}

double compute_mean(const std::vector<double> &values) {
    if (values.empty())
        return 0.0;
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    return sum / static_cast<double>(values.size());
}

double compute_variance(const std::vector<double> &values) {
    if (values.empty())
        return 0.0;
    double mean = compute_mean(values);
    double var_sum = 0.0;
    for (double x : values) {
        double diff = x - mean;
        var_sum += diff * diff;
    }
    return var_sum / static_cast<double>(values.size());
}
double compute_stddev(const std::vector<double> &values) { return std::sqrt(compute_variance(values)); }

std::vector<double> equally_spaced_points(int n, double start, double end) {
    std::vector<double> points;
    if (n <= 0)
        return points;

    if (n == 1) {
        points.push_back(start);
        return points;
    }

    points.reserve(n);
    double step = (end - start) / (n - 1); // n-1 intervals
    for (int i = 0; i < n; ++i) {
        points.push_back(start + i * step);
    }
    return points;
}

} // namespace math_utils

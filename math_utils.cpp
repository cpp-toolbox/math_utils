#include "math_utils.hpp"
#include <cmath>
#include <algorithm>

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

} // namespace math_utils

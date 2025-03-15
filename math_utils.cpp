#include "math_utils.hpp"
#include <cmath>

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

} // namespace math_utils

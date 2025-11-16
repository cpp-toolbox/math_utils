#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

#include <glm/glm.hpp>

namespace math_utils {

int non_neg_mod(int value, int mod);
std::pair<float, float> extract_yaw_pitch(const glm::vec3 &forward);

double map_range(double value, double in_min, double in_max, double out_min, double out_max);

} // namespace math_utils

#endif // MATH_UTILS_HPP

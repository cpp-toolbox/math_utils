#include "math_utils.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace math_utils {

int non_neg_mod(int value, int mod) { return (value % mod + mod) % mod; }

std::vector<glm::ivec2> get_square_grid_indices_along_line(int col1, int row1, int col2, int row2,
                                                           bool connect_line_segments_from_center_of_cell) {
    std::vector<glm::ivec2> grid_indices;

    if (!connect_line_segments_from_center_of_cell) {
        // Standard integer Bresenham
        int x = col1;
        int y = row1;

        int dx = std::abs(col2 - col1);
        int dy = std::abs(row2 - row1);

        int sx = (col1 < col2) ? 1 : -1;
        int sy = (row1 < row2) ? 1 : -1;

        int err = dx - dy;

        while (true) {
            grid_indices.push_back(glm::ivec2(x, y));

            if (x == col2 && y == row2)
                break;

            int e2 = 2 * err;

            if (e2 > -dy) {
                err -= dy;
                x += sx;
            }
            if (e2 < dx) {
                err += dx;
                y += sy;
            }
        }

    } else {
        // Floating-point DDA from cell centers
        float x0 = col1 + 0.5f;
        float y0 = row1 + 0.5f;
        float x1 = col2 + 0.5f;
        float y1 = row2 + 0.5f;

        float dx = x1 - x0;
        float dy = y1 - y0;

        int steps = static_cast<int>(std::max(std::abs(dx), std::abs(dy)));

        float x_inc = dx / steps;
        float y_inc = dy / steps;

        float x = x0;
        float y = y0;

        for (int i = 0; i <= steps; ++i) {
            grid_indices.push_back(glm::ivec2(static_cast<int>(x), static_cast<int>(y)));
            x += x_inc;
            y += y_inc;
        }
    }

    return grid_indices;
}

std::vector<glm::ivec2> get_square_grid_indices_in_annulus(int center_col, int center_row, float inner_radius,
                                                           float outer_radius) {
    std::vector<glm::ivec2> grid_indices;

    int min_col = static_cast<int>(std::floor(center_col - outer_radius));
    int max_col = static_cast<int>(std::ceil(center_col + outer_radius));
    int min_row = static_cast<int>(std::floor(center_row - outer_radius));
    int max_row = static_cast<int>(std::ceil(center_row + outer_radius));

    float center_x = center_col + 0.5f;
    float center_y = center_row + 0.5f;

    float inner_radius_sq = inner_radius * inner_radius;
    float outer_radius_sq = outer_radius * outer_radius;

    for (int row = min_row; row <= max_row; ++row) {
        for (int col = min_col; col <= max_col; ++col) {
            // center of current cell
            float cell_x = col + 0.5f;
            float cell_y = row + 0.5f;

            // distance from center cell to current cell (center-to-center)
            float dx = cell_x - center_x;
            float dy = cell_y - center_y;
            float dist_sq = dx * dx + dy * dy;

            // check if within annulus (inner_radius < distance <= outer_radius)
            if (dist_sq > inner_radius_sq && dist_sq <= outer_radius_sq) {
                grid_indices.push_back(glm::ivec2(col, row));
            }
        }
    }

    return grid_indices;
}

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

glm::vec3 hermite_spline(const glm::vec3 &p0, const glm::vec3 &p1, const glm::vec3 &t0, const glm::vec3 &t1, float t,
                         bool scale_tangents_by_distance) {
    float t2 = t * t;
    float t3 = t2 * t;

    float h00 = 2 * t3 - 3 * t2 + 1;
    float h10 = t3 - 2 * t2 + t;
    float h01 = -2 * t3 + 3 * t2;
    float h11 = t3 - t2;

    glm::vec3 scaled_t0 = t0;
    glm::vec3 scaled_t1 = t1;

    if (scale_tangents_by_distance) {
        float distance = glm::length(p1 - p0);
        scaled_t0 *= distance;
        scaled_t1 *= distance;
    }

    return h00 * p0 + h10 * scaled_t0 + h01 * p1 + h11 * scaled_t1;
}

glm::vec3 compute_best_fit_plane_normal(const std::vector<glm::vec3> &pts) {
    glm::vec3 normal(0.0f);

    const std::size_t count = pts.size();
    for (std::size_t i = 0; i < count; ++i) {
        const glm::vec3 &current = pts[i];
        const glm::vec3 &next = pts[(i + 1) % count];

        // NOTE: read about newell's method to understand what this is doing.
        normal.x += (current.y - next.y) * (current.z + next.z);
        normal.y += (current.z - next.z) * (current.x + next.x);
        normal.z += (current.x - next.x) * (current.y + next.y);
    }

    return normal;
}

bool points_are_planar(const std::vector<glm::vec3> &pts, float epsilon) {
    if (pts.size() < 3)
        return true;

    glm::vec3 normal = compute_best_fit_plane_normal(pts);
    float normal_len = glm::length(normal);

    if (normal_len < 1e-6f)
        return false; // degenerate / collinear

    normal /= normal_len;

    // use centroid as plane point
    glm::vec3 centroid(0.0f);
    for (const auto &p : pts)
        centroid += p;
    centroid /= static_cast<float>(pts.size());

    float max_distance = 0.0f;
    for (const auto &p : pts) {
        float distance = std::abs(glm::dot(p - centroid, normal));
        max_distance = std::max(max_distance, distance);
    }

    return max_distance <= epsilon;
}

} // namespace math_utils

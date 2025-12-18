#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

#include <algorithm>
#include <glm/glm.hpp>
#include <numeric>
#include <random>
#include <vector>
#include <deque>

namespace math_utils {

/**
 * @brief Computes the non-negative modulo of an integer.
 *
 * This function ensures that the result of modulo is always non-negative,
 * even if the input value is negative.
 *
 * @param value The integer value to reduce modulo.
 * @param mod The modulus.
 * @return The non-negative remainder of value divided by mod.
 */
int non_neg_mod(int value, int mod);

/**
 * @brief Extracts yaw and pitch angles (in degrees) from a 3D forward vector.
 *
 * The input vector is assumed to be a direction in 3D space. The function
 * normalizes the vector and computes:
 * - pitch: rotation around the X-axis (up/down)
 * - yaw: rotation around the Y-axis (left/right)
 *
 * @param forward A 3D forward vector (glm::vec3).
 * @return A pair of floats {yaw, pitch}, in degrees.
 *
 * @todo move to linalg utils
 */
std::pair<float, float> extract_yaw_pitch(const glm::vec3 &forward);

/**
 * @brief Maps a value from one numerical range to another.
 *
 * The input value is first clamped to [in_min, in_max], then linearly
 * interpolated to the output range [out_min, out_max].
 *
 * @param value The input value to map.
 * @param in_min Minimum of the input range.
 * @param in_max Maximum of the input range.
 * @param out_min Minimum of the output range.
 * @param out_max Maximum of the output range.
 * @return The value mapped to the output range.
 */
double map_range(double value, double in_min, double in_max, double out_min, double out_max);

/**
 * @brief Computes the arithmetic mean of a vector of doubles.
 *
 * If the input vector is empty, returns 0.0.
 *
 * @param values A vector of double values.
 * @return The mean of the values.
 */
double compute_mean(const std::vector<double> &values);

/**
 * @brief Computes the variance of a vector of doubles.
 *
 * Variance measures the average squared deviation from the mean.
 * If the input vector is empty, returns 0.0.
 *
 * @param values A vector of double values.
 * @return The variance of the values.
 */
double compute_variance(const std::vector<double> &values);

/**
 * @brief Computes the standard deviation of a vector of doubles.
 *
 * Standard deviation is the square root of the variance and provides
 * a measure of spread in the same units as the input values.
 *
 * @param values A vector of double values.
 * @return The standard deviation of the values.
 */
double compute_stddev(const std::vector<double> &values);

std::vector<double> equally_spaced_points(int n, double start, double end);

/*
 * @brief Exponential Moving Average (EMA) helper class.
 *
 * Tracks a moving average of incoming samples, weighting recent samples
 * more heavily and letting older samples decay exponentially.
 * Useful for tracking signals that vary over time without storing full history.
 *
 * The EMA is updated according to:
 * @code
 * EMA_new = alpha * new_sample + (1-alpha) * EMA_old
 * @endcode
 * This ensures that recent samples have more influence, and older samples' impact
 * diminishes exponentially over time.
 *
 * @param value The new sample value to include in the moving average.
 *
 * @example
 * MovingAverage avg(0.2); // alpha = 0.2 (so more smoothing)
 * avg.add_sample(10.0); // EMA = 10.0 (first sample)
 * avg.add_sample(20.0); // EMA = 0.2*20 + 0.8*10 = 12.0 (notice not a huge jump, aka smooth)
 * avg.add_sample(30.0); // EMA = 0.2*30 + 0.8*12 = 15.6
 * avg.add_sample(30.0); // EMA = 0.2*30 + 0.8*15.6 = 18.48
 * // Notice how the impact of the first sample (10.0) is diminishing over time
 */
class ExponentialMovingAverage {
  public:
    /// @brief Construct a MovingAverage with a given smoothing factor.
    /// @param alpha Smoothing factor in range (0,1]. Higher alpha = more responsive, lower alpha = smoother.
    explicit ExponentialMovingAverage(double alpha = 0.1) : alpha(alpha), has_value(false), current_average(0.0) {}

    /**
     * @brief Feed a new sample into the moving average.
     *
     */
    void add_sample(double value) {
        if (!has_value) {
            current_average = value;
            has_value = true;
        } else {
            current_average = (1.0 - alpha) * current_average + alpha * value;
        }
    }

    /**
     * @brief Get the current value of the moving average.
     * @return The current EMA value.
     */
    double get() const { return current_average; }

  private:
    double alpha;           ///< smoothing factor (0 < alpha <= 1)
    bool has_value;         ///< true if the first sample has been added
    double current_average; ///< current EMA value
};

class SimpleMovingAverage {
  public:
    /// @brief Construct a SMA with a fixed window size.
    /// @param window_size Number of recent samples to average.
    explicit SimpleMovingAverage(size_t window_size) : window_size(window_size), sum(0.0) {}

    /// @brief Add a new sample to the moving average.
    void add_sample(double value) {
        samples.push_back(value);
        sum += value;

        // Remove oldest sample if we exceed window size
        if (samples.size() > window_size) {
            sum -= samples.front();
            samples.pop_front();
        }
    }

    /// @brief Get the current simple moving average.
    double get() const {
        if (samples.empty())
            return 0.0;
        return sum / samples.size();
    }

  private:
    size_t window_size;         ///< number of samples to average
    std::deque<double> samples; ///< circular buffer for last N samples
    double sum;                 ///< sum of samples for fast average calculation
};

class PerlinNoise {
    // future: https://noiseposti.ng/posts/2022-01-16-The-Perlin-Problem-Moving-Past-Square-Noise.html
  public:
    explicit PerlinNoise(unsigned int seed = 2025) { init(seed); }

    std::vector<int> permutation_vector_for_hashing;

    void init(unsigned int seed) {
        permutation_vector_for_hashing.resize(256);
        std::iota(permutation_vector_for_hashing.begin(), permutation_vector_for_hashing.end(), 0);

        std::default_random_engine engine(seed);
        std::shuffle(permutation_vector_for_hashing.begin(), permutation_vector_for_hashing.end(), engine);

        // duplicate the permutation vector
        permutation_vector_for_hashing.insert(permutation_vector_for_hashing.end(),
                                              permutation_vector_for_hashing.begin(),
                                              permutation_vector_for_hashing.end());
    }

    /**
     * @brief Fade function used in Perlin noise to smooth interpolation.
     *
     * This function is a quintic polynomial that eases values in the range [0, 1].
     * It ensures smooth transitions between lattice points in Perlin noise by
     * providing zero first and second derivatives at the endpoints.
     *
     * The function is defined as:
     * @f[
     *   \text{fade}(t) = 6t^5 - 15t^4 + 10t^3
     * @f]
     *
     * Properties:
     * - fade(0) = 0
     * - fade(1) = 1
     * - fade'(0) = fade'(1) = 0
     * - fade''(0) = fade''(1) = 0
     *
     * These properties produce a smooth "S-curve" that eases in and out,
     * making linear interpolation between noise values visually smooth and
     * continuous.
     *
     * @param t A value in the range [0, 1] representing the interpolation factor.
     * @return The eased value corresponding to t.
     */
    inline double fade(double t) const { return t * t * t * (t * (t * 6 - 15) + 10); }

    // linear interpolation
    inline double lerp(double t, double a, double b) const { return a + t * (b - a); }

    /**
     * @brief Gradient function used in Perlin noise to calculate directional vectors.
     *
     * This function determines a pseudo-random gradient vector based on a hash value
     * and computes the dot product of that vector with the distance vector (x, y, z)
     * from the lattice point. It is a key step in generating smooth Perlin noise.
     *
     * How it works:
     * - The hash value is masked with 15 to keep only the lower 4 bits: `h = hash & 15`.
     * - This 4-bit value selects one of 16 possible gradient directions in 3D space.
     * - `u` and `v` are chosen components of (x, y, z) based on `h`:
     *   - `u = x` if h < 8, otherwise `u = y`
     *   - `v = y` if h < 4, `v = x` if h == 12 or 14, otherwise `v = z`
     * - The function returns a combination of `u` and `v` with signs determined by the bits in `h`:
     *   - If bit 0 of h is set, `u` is negated
     *   - If bit 1 of h is set, `v` is negated
     * - Resulting value is equivalent to the dot product of the chosen gradient vector and the input vector.
     *
     * This design avoids storing explicit gradient vectors while producing pseudo-random
     * gradients that are consistent and repeatable, which is critical for smooth noise.
     *
     * @param hash An integer hash value, usually derived from a permutation table.
     * @param x The x-component of the distance vector from the lattice point.
     * @param y The y-component of the distance vector from the lattice point.
     * @param z The z-component of the distance vector from the lattice point.
     * @return The dot product of the gradient vector selected by the hash and the input vector (x, y, z).
     */
    inline double grad(int hash, double x, double y, double z) const {
        int h = hash & 15; // Take the lower 4 bits of the hash
        double u = h < 8 ? x : y;
        double v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
        return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
    }

    double at(double x, double y, double z) const {
        // integer parts (wrapped to 0-255)
        int floor_x = static_cast<int>(floor(x)) & 255;
        int floor_y = static_cast<int>(floor(y)) & 255;
        int floor_z = static_cast<int>(floor(z)) & 255;

        double frac_x = x - floor(x);
        double frac_y = y - floor(y);
        double frac_z = z - floor(z);

        // fade curves
        double u = fade(frac_x);
        double v = fade(frac_y);
        double w = fade(frac_z);

        // hash coordinates of the cube corners
        int hash_XY = permutation_vector_for_hashing[floor_x] + floor_y;
        int hash_X1Y = permutation_vector_for_hashing[floor_x + 1] + floor_y;

        int corner000 = permutation_vector_for_hashing[hash_XY] + floor_z;         // AA
        int corner001 = permutation_vector_for_hashing[hash_XY] + floor_z + 1;     // AA + 1
        int corner010 = permutation_vector_for_hashing[hash_XY + 1] + floor_z;     // AB
        int corner011 = permutation_vector_for_hashing[hash_XY + 1] + floor_z + 1; // AB + 1

        int corner100 = permutation_vector_for_hashing[hash_X1Y] + floor_z;         // BA
        int corner101 = permutation_vector_for_hashing[hash_X1Y] + floor_z + 1;     // BA + 1
        int corner110 = permutation_vector_for_hashing[hash_X1Y + 1] + floor_z;     // BB
        int corner111 = permutation_vector_for_hashing[hash_X1Y + 1] + floor_z + 1; // BB + 1

        double grad000 = grad(permutation_vector_for_hashing[corner000], frac_x, frac_y, frac_z);
        double grad100 = grad(permutation_vector_for_hashing[corner100], frac_x - 1, frac_y, frac_z);
        double grad010 = grad(permutation_vector_for_hashing[corner010], frac_x, frac_y - 1, frac_z);
        double grad110 = grad(permutation_vector_for_hashing[corner110], frac_x - 1, frac_y - 1, frac_z);

        double grad001 = grad(permutation_vector_for_hashing[corner001], frac_x, frac_y, frac_z - 1);
        double grad101 = grad(permutation_vector_for_hashing[corner101], frac_x - 1, frac_y, frac_z - 1);
        double grad011 = grad(permutation_vector_for_hashing[corner011], frac_x, frac_y - 1, frac_z - 1);
        double grad111 = grad(permutation_vector_for_hashing[corner111], frac_x - 1, frac_y - 1, frac_z - 1);

        // interpolate along x
        double lerp_x00 = lerp(u, grad000, grad100);
        double lerp_x10 = lerp(u, grad010, grad110);
        double lerp_x01 = lerp(u, grad001, grad101);
        double lerp_x11 = lerp(u, grad011, grad111);

        // interpolate along y
        double lerp_y0 = lerp(v, lerp_x00, lerp_x10);
        double lerp_y1 = lerp(v, lerp_x01, lerp_x11);

        // interpolate along z
        double result = lerp(w, lerp_y0, lerp_y1);

        // normalize to [0,1]
        return (result + 1.0) / 2.0;
    }

    /**
     * @brief this function returns an image of the current noise
     *
     * usage:
     * @code
     * size_t width = 512, height = 512;
     * auto image = perlin_noise.get_noise_image(width, height);
     * stbi_write_png("perlin.png", width, height, 1, image.data(), width);
     * @endcode
     *
     */
    std::vector<unsigned char> get_noise_image(size_t width = 512, size_t height = 512) {

        std::vector<unsigned char> image(width * height);

        double scale = 0.01; // zoom in/out of noise

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                double n = at(x * scale, y * scale, 0.0); // 0..1
                image[y * width + x] = static_cast<unsigned char>(n * 255);
            }
        }

        return image;
    }
};

} // namespace math_utils

#endif // MATH_UTILS_HPP

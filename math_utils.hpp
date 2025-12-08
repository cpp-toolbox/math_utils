#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP

#include <glm/glm.hpp>
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

} // namespace math_utils

#endif // MATH_UTILS_HPP

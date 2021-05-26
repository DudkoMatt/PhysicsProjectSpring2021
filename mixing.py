from main import *
import scipy.optimize


def mix_it_2(first_color: list[float], second_color: list[float], first_color_part, second_color_part) -> list[float]:
    """
    :param first_color: first color to be mixed
    :param second_color: second color to be mixed
    color list represents normalised from 0 to 1 intensity of wavelength
    :param first_color_part: part of first color
    :param second_color_part: part of second color
    Part are numbers, so that mixing will be in first_color_part : second_color_part proportion
    :return: list of intensity normalised from 0 to 1 of mixed color
    """
    _all = first_color_part + second_color_part
    first = first_color_part / _all
    second = second_color_part / _all
    result = []
    assert len(first_color) == len(second_color)
    for i in range(0, len(first_color)):
        mixed = (first_color[i] ** first) * (second_color[i] ** second)
        result.append(mixed)
    return result


def mix_it(colors: list[list[float]], colors_parts: list[float]) -> list[float]:
    """
    :param colors: list of colors to mix
    color list represents normalised from 0 to 1 intensity of wavelength
    :param colors_parts: list of parts of each color
    For example colors = ['red', 'green', 'blue'], colors_parts = [1, 2, 3]
    In the result will be 1/6     of red
                          2/6=1/3 of green
                          3/6=1/2 of blue
    :return: list of intensity normalised from 0 to 1 of mixed color
    """
    total = sum(colors_parts)
    part_coefficients = list(map(lambda _x: _x / total, colors_parts))
    result = [1] * len(colors[0])

    # Check all lengths of colors
    for color_idx in range(len(colors) - 1):
        assert len(colors[color_idx]) == len(colors[color_idx + 1]), "Lengths of specters does not match"

    for _i in range(0, len(colors)):
        for _j in range(0, len(colors[_i])):
            result[_j] *= colors[_i][_j] ** part_coefficients[_i]
    return result


def spectrum_diff(first_spectrum: list[float], second_spectrum: list[float]) -> float:
    assert(len(first_spectrum) == len(second_spectrum))

    return math.sqrt(sum([(first_spectrum[_i] - second_spectrum[_i]) ** 2 for _i in range(len(first_spectrum))]))


def analyze_color(base_colors_spectrum: list[list[float]], color_to_analyze: list[float]) -> list[float]:
    """
    Main function to analyze unknown color
    :param base_colors_spectrum: spectrum of base colors to mix
    :param color_to_analyze: spectrum of color to analyze
    :return: proportion coefficients to mix base colors
    """

    colors_parts_coefficients = np.array([1] * len(base_colors_spectrum))

    def function_to_optimize(k_coefficients: np.array):
        current_calculated_color = mix_it(base_colors_spectrum, k_coefficients)
        return spectrum_diff(current_calculated_color, color_to_analyze)

    result = scipy.optimize.minimize(function_to_optimize, colors_parts_coefficients)
    minimum = min(result.x)
    return list(map(lambda element: round(element / minimum, 5), result.x))

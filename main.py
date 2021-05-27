import os
from scipy.signal import savgol_filter
from scipy import interpolate as interp
import matplotlib.pyplot
import numpy as np
from numpy import linalg as LA
import re
from os import listdir, getcwd
from os.path import isfile, join
from typing import Callable
import sys
import platform
from PIL import Image
import math


# Getting vectors from file
def get_all_vectors_from_file(filename: str, naming_pattern=r"^\[[0-9]*\]-.*-[0-9]*_E$") -> [list[list[float]],
                                                                                             list[float]]:
    """
    Filtering data from *.lab
    :param filename: name of the file
    :param naming_pattern: which vectors not to ignore
    :return: all vectors, lambda
    """

    with open(filename, 'r', encoding='utf16') as f:
        total_str = ''.join(f.readlines())

    lambda_was_already = False
    lambda_array = []
    vectors = []

    while total_str.find('vecteur') != -1:
        # Skip string to next vector
        total_str = total_str[total_str.find('vecteur') + len('vecteur]\n\t'):]
        total_str = total_str[total_str.find('\n\t') + 2:]

        # Get name
        name = total_str[:total_str.find('\n')]
        name = name[name.find('"') + 1:-1]

        # Getting data
        total_str = total_str[total_str.find('points'):]
        data = total_str[:total_str.find("}") + 1]
        data = data[data.find("{") + 1:-1].replace("\t", " ").replace("\n", " ").split()

        # Ignore RAW and lambda (if more than 1 time)

        if name == 'λ' and lambda_was_already:
            continue

        pattern = re.compile(naming_pattern)

        if not pattern.match(name) and name != 'λ':
            # Ignoring
            continue

        if name == 'λ':
            lambda_was_already = True
            lambda_array = list(map(float, data))
            continue

        vectors.append(list(map(float, data)))

    start = np.where(np.array(lambda_array) >= 410)[0][0]
    end_gt_410 = np.where(np.array(lambda_array) >= 410)[0][-1] + 1
    end = np.where(np.array(lambda_array[:end_gt_410]) <= 700)[0][-1] + 1

    return [vector[start:end] for vector in vectors], lambda_array[start:end]


def average(arrays: list[list[float]]) -> list[float]:  # ToDO: refactor
    """
    Get the average for each row in 2-d array
    :param arrays: 2-d array, each list (column) is one measure
    :return: 1-d array with average for each row
    """

    result = []
    for i in range(len(arrays[0])):
        count = len(arrays)
        sum = 0
        for array in arrays:
            sum += array[i]
        av = sum / count
        result.append(av)
    return result


def normalize(vector: list[float], white_vector: list[float], black_vector: list[float]) -> list[float]:
    """
    Normalize vector data to range [0, 1] relative to white vector data and black vector data
    :param vector: vector data (one measurement)
    :param white_vector: vector data of white color
    :param black_vector: vector data of black color
    :return: normalized vector data relative to white and black colors
    """

    result = []
    for idx in range(len(vector)):
        result.append((vector[idx] - black_vector[idx]) / (white_vector[idx] - black_vector[idx]))

    return result


def median(array: list[float]) -> list[float]:
    """
    Apply median filtering and cut vector data (reduced by 4 elements)
    :param array: vector data (one measurement)
    :return: vector data after filtering
    """

    result = []
    for i in range(2, len(array) - 2):
        result.append(np.median(array[i - 2:i + 3]))

    return result


def interpolate(array: list[float], window_length=151, polyorder=1) -> np.ndarray:
    """
    Apply filtering and interpolating (savgol_filter) to vector data
    :param array: vector data (one measurement)
    :param window_length: parameter for savgol_filter
    :param polyorder: parameter for savgol_filter
    :return: vector data after filtering and interpolating
    """

    return savgol_filter(array, window_length, polyorder)


def plot(x: list[float], y: np.ndarray, plot_name: str, normalised=True):
    """
    Plotting data with matplotlib
    :param x: x values vector data
    :param y: y values vector data
    :param plot_name: the name of the plot
    :param normalised: is data normalized to [0, 1] or not. Default is True
    :return: None
    """

    return   # -----------------------------------------------------TODO: REMOVE!----------------------------------------------------------------------------------------------------
    matplotlib.pyplot.plot(x, y)
    matplotlib.pyplot.title(plot_name)
    matplotlib.pyplot.xlabel('λ, nm')
    matplotlib.pyplot.ylabel('E_0' if normalised else 'E')
    if normalised:
        matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.show()


def analyze(all_data: list[list[list[float]]], lambda_data: list[float], data_namings: list[str], white_idx: int,
            black_idx: int) -> list[list[float]]:
    """
    Main function for analyzing, filtering, interpolating, plotting data
    :param all_data: several samples data 3-d array (different samples -> several measurements -> measurement)
    :param lambda_data: lambda values for which measurements are taken
    :param data_namings: names for the samples
    :param white_idx: white data index in all_data list
    :param black_idx: black data index in all_data list
    :return: data after calculations
    """

    white_data_average = average(all_data[white_idx])
    black_data_average = average(all_data[black_idx])

    result = []

    for idx in range(len(all_data)):
        # one_file_data = interpolate(all_data[idx])
        one_file_data = average(all_data[idx])
        # one_file_data = average(one_file_data)

        one_file_data = median(one_file_data)

        if not (idx == white_idx or idx == black_idx):
            one_file_data = normalize(one_file_data, white_data_average, black_data_average)  # ToDO: interpolated?

        one_file_data = interpolate(one_file_data)

        plot_name = data_namings[idx]
        normalized = True
        if idx == white_idx:
            plot_name = 'White - ' + plot_name
            normalized = False
        elif idx == black_idx:
            plot_name = 'Black - ' + plot_name
            normalized = False

        result.append(one_file_data)
        plot(lambda_data, one_file_data, plot_name, normalized)

    return result


class InterpolatedColorCoefficients:
    """
    Class for calculation color coefficients used for getting values from spectrometer data

    Attributes
    ----------
    filename: str
        a name of the file from which retrieve data
    color_coefficients: dict
        data from file with color coefficients
    x_func: Callable[[float], float]
        interpolation function for x part
    y_func: Callable[[float], float]
        interpolation function for y part
    z_func: Callable[[float], float]
        interpolation function for z part
    """

    filename: str
    color_coefficients: dict  # lambda -> (x, y, z)
    x_func: Callable[[float], float]
    y_func: Callable[[float], float]
    z_func: Callable[[float], float]

    def __init__(self, filename='coefficients.txt') -> None:
        self.filename = filename
        self.color_coefficients = self.get_color_coefficients(filename)
        self.x_func, self.y_func, self.z_func = self.__interpolate_color_coefficients()
        # pass

    @staticmethod
    def get_color_coefficients(filename: str) -> dict:
        # lambda -> (x, y, z)
        with open(filename, 'r') as f:
            coefficients = list(
                map(lambda item: list(map(float, item.replace('\n', '').replace(',', '.').split())), f.readlines()))

        result_dict = dict()
        for _lambda, x, y, z in coefficients:
            result_dict[_lambda] = (x, y, z)

        return result_dict

    def __interpolate_color_coefficients(self) -> (
    Callable[[float], float], Callable[[float], float], Callable[[float], float]):
        lambda_array = []
        x_array = []
        y_array = []
        z_array = []

        for key in self.color_coefficients.keys():
            lambda_array.append(key)
            x_array.append(self.color_coefficients[key][0])
            y_array.append(self.color_coefficients[key][1])
            z_array.append(self.color_coefficients[key][2])

        return \
            interp.interp1d(lambda_array, x_array), \
            interp.interp1d(lambda_array, y_array), \
            interp.interp1d(lambda_array, z_array)

    def get_coefficients_for(self, wavelength: float) -> (float, float, float):
        return Coefficient(self.x_func(wavelength), self.y_func(wavelength), self.z_func(wavelength))


class Coefficient:
    x: float
    y: float
    z: float

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def __mul__(self, other):
        assert isinstance(other, float)

        self.x *= other
        self.y *= other
        self.z *= other

        return self

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        assert isinstance(other, Coefficient)

        self.x += other.x
        self.y += other.y
        self.z += other.z

        return self

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        return self.__add__(other)

    def __imul__(self, other):
        return self.__add__(other)


def sum_of_one_vector_data(data: list[float], lambda_data: list[float],
                           color_coefficients: InterpolatedColorCoefficients, DELTA_LAMBDA: float) -> Coefficient:
    result = Coefficient(0, 0, 0)
    length = len(lambda_data)
    for i in range(0, length):
        result += data[i] * color_coefficients.get_coefficients_for(lambda_data[i]) * DELTA_LAMBDA
    return result


def calculate_color_xyz(data_vector: list[float], lambda_data: list[float],
                        color_coefficients: InterpolatedColorCoefficients, DELTA_LAMBDA: float) -> (float, float, float):
    coefficient = sum_of_one_vector_data(data_vector, lambda_data, color_coefficients, DELTA_LAMBDA)
    k_c = 100 / coefficient.y

    x = k_c * coefficient.x
    y = 100
    z = k_c * coefficient.z

    return (X := x / (x + y + z)), (Y := y / (x + y + z)), 1 - X - Y


class TransitionCoefficient:
    x: float
    y: float
    z: float

    def __init__(self, x: float, y: float) -> None:
        self.x = x / y
        self.y = 1
        self.z = (1 - x - y) / y


def get_A_white_coefficient():
    x_w = 109.850 / 100  # коэффициенты взяты со стр 110 для источника A
    y_w = 100 / 100
    z_w = 35.585 / 100
    return np.array([[x_w], [y_w], [z_w]])


def get_E_white_coefficient():
    x_w = 100 / 100  # коэффициенты взяты со стр 110 для источника E
    y_w = 100 / 100
    z_w = 100 / 100
    return np.array([[x_w], [y_w], [z_w]])


def get_D65_white_coefficient():
    x_w = 95.047 / 100  # коэффициенты взяты со стр 110 для источника D65
    y_w = 100 / 100
    z_w = 108.883 / 100
    return np.array([[x_w], [y_w], [z_w]])


def get_M_A_matrix(xyz_scaling=False, bradford=False, von_kries=False):
    """
        # XYZ Scaling
        # -----------
        M_A = np.array([[1.0000000 0.0000000 0.0000000],
        [0.0000000 1.0000000 0.0000000],
        [0.0000000 0.0000000 1.0000000]])

        M_A_inv = np.array([[1.0000000 0.0000000 0.0000000],
        [0.0000000 1.0000000 0.0000000],
        [0.0000000 0.0000000 1.0000000]])


        # Bradford
        # --------
        M_A = np.array([[0.8951000 0.2664000 -0.1614000],
        [-0.7502000 1.7135000 0.0367000],
        [0.0389000 -0.0685000 1.0296000]])

        M_A_inv = np.array([[0.9869929 -0.1470543 0.1599627],
        [0.4323053 0.5183603 0.0492912],
        [-0.0085287 0.0400428 0.9684867]])


        # Von Kries
        # ---------
        M_A = np.array([[0.4002400 0.7076000 -0.0808100],
        [-0.2263000 1.1653200 0.0457000],
        [0.0000000 0.0000000 0.9182200]])

        M_A_inv = np.array([[1.8599364 -1.1293816 0.2198974],
        [0.3611914 0.6388125 -0.0000064],
        [0.0000000 0.0000000 1.0890636]])
    """
    if xyz_scaling:
        return np.array([[1.0000000, 0.0000000, 0.0000000],
                    [0.0000000, 1.0000000, 0.0000000],
                    [0.0000000, 0.0000000, 1.0000000]])
    elif bradford:
        return np.array([[0.8951000, 0.2664000, -0.1614000],
                        [-0.7502000, 1.7135000, 0.0367000],
                        [0.0389000, -0.0685000, 1.0296000]])
    elif von_kries:
        return np.array([[0.4002400, 0.7076000, -0.0808100],
                        [-0.2263000, 1.1653200, 0.0457000],
                        [0.0000000, 0.0000000, 0.9182200]])
    else:
        print("Wow, you should pass an argument!", file=sys.stderr)
        exit(-1)


def convert_white_point(x: float, y: float, z: float, from_white: np.array, to_white: np.array) -> (float, float, float):
    M_A = get_M_A_matrix(von_kries=True)
    M_A_inv = LA.inv(M_A)

    from_correction_coefs = M_A.dot(from_white)
    to_correction_coefs = M_A.dot(to_white)

    first = to_correction_coefs[0, 0] / from_correction_coefs[0, 0]
    second = to_correction_coefs[1, 0] / from_correction_coefs[1, 0]
    third = to_correction_coefs[2, 0] / from_correction_coefs[2, 0]

    M_C = M_A_inv.dot(np.diag([first, second, third])).dot(M_A)

    return M_C.dot(np.array([[x], [y], [z]]))


def get_inverse_m_matrix():
    x_r = 0.7350  # Коэффициенты из приложения 4 книжки на стр 111 тип CIE RGB
    y_r = 0.2650
    x_g = 0.2740
    y_g = 0.7170
    x_b = 0.1670
    y_b = 0.0090

    r_avg = TransitionCoefficient(x_r, y_r)
    g_avg = TransitionCoefficient(x_g, y_g)
    b_avg = TransitionCoefficient(x_b, y_b)
    s = (LA.inv(np.array([[r_avg.x, g_avg.x, b_avg.x], [r_avg.y, g_avg.y, b_avg.y], [r_avg.z, g_avg.z, b_avg.z]]))).dot(
        get_E_white_coefficient())
    s_r = s[0, 0]
    s_g = s[1, 0]
    s_b = s[2, 0]

    return LA.inv(np.array(
        [[s_r * r_avg.x, s_g * g_avg.x, s_b * b_avg.x],
         [s_r * r_avg.y, s_g * g_avg.y, s_b * b_avg.y],
         [s_r * r_avg.z, s_g * g_avg.z, s_b * b_avg.z]]))


def xyz_to_linear_rgb(x: float, y: float, z: float) -> (float, float, float):
    result_vector = get_inverse_m_matrix().dot(np.array([x, y, z]))
    return result_vector[0], result_vector[1], result_vector[2]


def rgb_color_diff(rgb1: np.array, rgb2: np.array) -> float:
    assert rgb1.shape == (3,), "Data should be in one row, not matrix 3x1"
    assert rgb2.shape == (3,), "Data should be in one row, not matrix 3x1"
    assert all([isinstance(_i, np.int32) and 0 <= _i <= 255 for _i in rgb1]), "Data should be in range 0..255, not 0..1 or negative"
    assert all([isinstance(_i, np.int32) and 0 <= _i <= 255 for _i in rgb2]), "Data should be in range 0..255, not 0..1 or negative"

    return math.sqrt(sum([(rgb1[_i] - rgb2[_i]) ** 2 for _i in range(len(rgb1))]))


# CAUTION: USING GLOBAL VARIABLES
def spectrum_to_rgb(spectrum: list[float]):
    _xyz_coord = calculate_color_xyz(spectrum, lambda_data, color_coefficients, DELTA_LAMBDA)
    return "{} {} {}".format(*list(map(lambda _x: max(min(round(_x * 255), 255), 0), xyz_to_linear_rgb(*_xyz_coord))))


if __name__ == "__main__":
    assert sys.version_info >= (3, 9)

    # Const from .lab file
    DELTA_LAMBDA = 0.086

    # Black and white indexes
    WHITE_IDX = 0
    BLACK_IDX = 1

    color_coefficients = InterpolatedColorCoefficients('coefficients.txt')

    # Getting files *.lab in current working directory

    os.chdir("Data_from_12.04.2021")
    os.chdir("Try #3")

    only_files = [f for f in listdir(getcwd()) if isfile(join(getcwd(), f))]
    files = list(filter(lambda x: re.match(r'.*\.lab$', x), only_files))

    all_data = []
    filenames = []

    naming_pattern = re.compile(r'([0-9]+)\.lab$')
    lambda_data = []

    for i in sorted(files, key=lambda x: naming_pattern.search(x).group().replace('.lab', '')):
        data_array, lambda_data = get_all_vectors_from_file(i)
        all_data.append(data_array)
        filenames.append(naming_pattern.search(i).group().replace('.lab', ''))

    # Cut lambda array due to using median
    lambda_data = lambda_data[2: len(lambda_data) - 2]

    processed_data = analyze(all_data, lambda_data, filenames, WHITE_IDX, BLACK_IDX)
    calculated_colors = []

    if platform.system() == "Windows":
        xyz_file = open('xyz.txt', 'w')
        open('xyz_converted.txt', 'w').close()
    else:
        xyz_file = open('/Users/igorklyuzev/ITMO/4_семестр/Физика/PhysicsProjectSpring2021/xyz.txt', 'w')

    k = 0

    for vector in processed_data:
        xyz_coord = calculate_color_xyz(vector, lambda_data, color_coefficients, DELTA_LAMBDA)
        print("{:.5f} {:.5f} {:.5f}".format(*xyz_coord), filenames[k], file=xyz_file)
        k += 1
        
        calculated_colors.append(xyz_to_linear_rgb(*xyz_coord))
    if platform.system() == "Windows":
        output_file = open('output.txt', 'w')
        for calculated_color in calculated_colors:
            print("{} {} {}".format(*list(map(lambda _x: max(min(round(_x * 255), 255), 0), calculated_color))), file=output_file)

        output_file.flush()
        output_file.close()
    else:
        print(*calculated_colors, sep='\n', file=open('/Users/igorklyuzev/ITMO/4_семестр/Физика/PhysicsProjectSpring2021/output.txt', 'w'))  # ToDO: wrong answers -> negative rgb coordinates

    k = 0
    for color in calculated_colors:
        img = Image.new('RGB', (1000, 1000), (max(min(round(color[0] * 255), 255), 0), max(min(round(color[1] * 255), 255), 0), max(min(round(color[2] * 255), 255), 0)))
        img.save("./images/{}.png".format(filenames[k]))
        k += 1

    xyz_file.flush()
    xyz_file.close()

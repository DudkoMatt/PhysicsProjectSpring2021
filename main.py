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
    result = []
    for idx in range(len(vector)):
        result.append((vector[idx] - black_vector[idx]) / (white_vector[idx] - black_vector[idx]))

    return result


def median(array: list[float]) -> list[float]:
    result = []
    for i in range(2, len(array) - 2):
        result.append(np.median(array[i - 2:i + 3]))

    return result


def interpolate(array: list[float], window_length=151, polyorder=1) -> np.ndarray:
    return savgol_filter(array, window_length, polyorder)


def plot(x: list[float], y: np.ndarray, plot_name: str, normalised=True):
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
        # TODO: SWITCH ------------------------------------------------------------------------------------------------------------------------------------------
        return Coefficient(self.x_func(wavelength), self.y_func(wavelength), self.z_func(wavelength))
        # return Coefficient(*InterpolatedColorCoefficients.get_color_coefficients_for(wavelength))

    @staticmethod
    def gaussian(x, alpha, mu, sigma_1, sigma_2):
        return alpha * np.exp(-((x - mu) / (sigma_1 if x < mu else sigma_2)) ** 2 / 2)

    @staticmethod
    def x_coefficient_CIE_1931(wavelength):
        return InterpolatedColorCoefficients.gaussian(wavelength, 1.056, 599.8, 379, 310) \
               + InterpolatedColorCoefficients.gaussian(wavelength, 0.362, 442.0, 160, 267) \
               + InterpolatedColorCoefficients.gaussian(wavelength, -0.065, 501.1, 204, 262)

    @staticmethod
    def y_coefficient_CIE_1931(wavelength):
        return InterpolatedColorCoefficients.gaussian(wavelength, 0.821, 568.8, 469, 405) \
               + InterpolatedColorCoefficients.gaussian(wavelength, 0.286, 530.9, 163, 311)

    @staticmethod
    def z_coefficient_CIE_1931(wavelength):
        return InterpolatedColorCoefficients.gaussian(wavelength, 1.217, 437.0, 118, 360) \
               + InterpolatedColorCoefficients.gaussian(wavelength, 0.681, 459.0, 260, 138)

    @staticmethod
    def get_color_coefficients_for(wavelength):
        return InterpolatedColorCoefficients.x_coefficient_CIE_1931(wavelength), \
               InterpolatedColorCoefficients.y_coefficient_CIE_1931(wavelength), \
               InterpolatedColorCoefficients.z_coefficient_CIE_1931(wavelength)


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
                           color_coefficients: InterpolatedColorCoefficients) -> Coefficient:
    result = Coefficient(0, 0, 0)
    length = len(lambda_data)
    for i in range(0, length):
        result += data[i] * color_coefficients.get_coefficients_for(lambda_data[i]) * DELTA_LAMBDA
    return result


def calculate_color_xyz(data_vector: list[float], lambda_data: list[float],
                        color_coefficients: InterpolatedColorCoefficients) -> (float, float, float):
    coefficient = sum_of_one_vector_data(data_vector, lambda_data, color_coefficients)
    k_c = 100 / coefficient.y

    x = k_c * coefficient.x
    y = 100
    z = k_c * coefficient.z

    return x, y, z


class TransitionCoefficient:
    x: float
    y: float
    z: float

    def __init__(self, x: float, y: float) -> None:
        self.x = x / y
        self.y = 1
        self.z = (1 - x - y) / y


def get_white_coefficient():
    x_w = 109.850  # ToDo коэффициенты взяты со стр 110 для источника A
    y_w = 100
    z_w = 35.585
    return np.array([[x_w], [y_w], [z_w]])


def get_inverse_m_matrix():
    # x_r = 0.7350  # ToDo коэффициенты из приложения 4 книжки на стр 111 тип CIE RGB
    # y_r = 0.2650  # ToDo
    # x_g = 0.2740  # ToDo
    # y_g = 0.7170  # ToDo
    # x_b = 0.1670  # ToDo
    # y_b = 0.0090  # ToDo

    x_r = 0.6400  # ToDo коэффициенты из приложения 4 книжки на стр 111 тип Adobe RGB
    y_r = 0.3300  # ToDo
    x_g = 0.2100  # ToDo
    y_g = 0.7100  # ToDo
    x_b = 0.1500  # ToDo
    y_b = 0.0600  # ToDo
    r_avg = TransitionCoefficient(x_r, y_r)
    g_avg = TransitionCoefficient(x_g, y_g)
    b_avg = TransitionCoefficient(x_b, y_b)
    s = (LA.inv(np.array([[r_avg.x, g_avg.x, b_avg.x], [r_avg.y, g_avg.y, b_avg.y], [r_avg.z, g_avg.z, b_avg.z]]))).dot(
        get_white_coefficient())
    s_r = s[0, 0]
    s_g = s[1, 0]
    s_b = s[2, 0]
    return LA.inv(np.array(
        [[s_r * r_avg.x, s_g * g_avg.x, s_b * b_avg.x], [s_r * r_avg.y, s_g * g_avg.y, s_b * b_avg.y],
         [s_r * r_avg.z, s_g * g_avg.z, s_b * b_avg.z]]))


def xyz_to_linear_rgb(x: float, y: float, z: float) -> (float, float, float):
    return get_inverse_m_matrix().dot(np.array([[x], [y], [z]]))[0, 0], \
           get_inverse_m_matrix().dot(np.array([[x], [y], [z]]))[1, 0], \
           get_inverse_m_matrix().dot(np.array([[x], [y], [z]]))[2, 0]


def linear_rgb_to_non_linear_rgb(x: float, y: float, z: float) -> (float, float, float):
    x_linear, y_linear, z_linear = xyz_to_linear_rgb(x, y, z)
    gamma = 2.2  # ToDo
    return x_linear ** (1 / gamma), y_linear ** (1 / gamma), z_linear ** (1 / gamma)


if __name__ == "__main__":
    assert sys.version_info >= (3, 9)

    # Const from .lab file
    DELTA_LAMBDA = 0.086
    # ToDO: больше нет необхоимости
    # LAMBDA_START = 400.022
    # LAMBDA_STOP = 699.978

    # Black and white indexes
    WHITE_IDX = 0
    BLACK_IDX = 1

    color_coefficients = InterpolatedColorCoefficients('coefficients.txt')

    # Getting files *.lab in current working directory

    os.chdir("Data_from_12.04.2021\\Try #3\\")

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

    # ------------------------------ Тестирование интерполяции --------------------------------------------------------
    # # ToDO
    #
    # coefficients_list = [[], [], [], []]
    #
    # for key in color_coefficients.keys():
    #     coefficients_list[0].append(key)
    #     coefficients_list[1].append(color_coefficients[key][0])
    #     coefficients_list[2].append(color_coefficients[key][1])
    #     coefficients_list[3].append(color_coefficients[key][2])
    #
    # interpolated_x = np.interp(lambda_data, coefficients_list[0], coefficients_list[1])
    # interpolated_y = np.interp(lambda_data, coefficients_list[0], coefficients_list[2])
    # interpolated_z = np.interp(lambda_data, coefficients_list[0], coefficients_list[3])
    #
    # # ToDO: построить графики -> коэффициенты до интеполяции и после интерполяции (3 * 2 графика для каждой из координат)
    #
    # exit(0)  # ToDO: убрать
    # -----------------------------------------------------------------------------------------------------------------

    processed_data = analyze(all_data, lambda_data, filenames, WHITE_IDX, BLACK_IDX)
    calculated_colors = []

    xyz_file = open('xyz.txt', 'w')

    for vector in processed_data:
        xyz_coord = calculate_color_xyz(vector, lambda_data, color_coefficients)
        print(xyz_coord, file=xyz_file)
        calculated_colors.append(xyz_to_linear_rgb(*xyz_coord))

    print(*calculated_colors, sep='\n', file=open('output.txt', 'w'))  # ToDO: wrong answers -> negative rgb coordinates

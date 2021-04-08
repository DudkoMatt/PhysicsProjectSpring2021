from scipy.signal import savgol_filter
import matplotlib.pyplot
import numpy as np
import re
from os import listdir, getcwd
from os.path import isfile, join


# Getting vectors from file
def get_all_vectors_from_file(filename: str, naming_pattern=r"^\[[0-9]*\]-.*-[0-9]*_E$") -> [list[list[float]], list[float]]:
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

    start = np.where(np.array(lambda_array) >= 400)[0][0]
    end_gt_400 = np.where(np.array(lambda_array) >= 400)[0][-1] + 1
    end = np.where(np.array(lambda_array[:end_gt_400]) <= 700)[0][-1] + 1

    return [vector[start:end] for vector in vectors], lambda_array[start:end]


def average(arrays: list[list[float]]) -> list[float]:
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
        result.append(np.median(array[i-2:i+3]))

    return result


def interpolate(array: list[float], window_length=151, polyorder=1) -> np.ndarray:
    return savgol_filter(array, window_length, polyorder)


def plot(x: list[float], y: np.ndarray, plot_name: str, normalised=True):
    matplotlib.pyplot.plot(x, y)
    matplotlib.pyplot.title(plot_name)
    matplotlib.pyplot.xlabel('λ, nm')
    matplotlib.pyplot.ylabel('E_0' if normalised else 'E')
    if normalised:
        matplotlib.pyplot.ylim(0, 1)
    matplotlib.pyplot.show()


def analyze(all_data: list[list[list[float]]], lambda_data: list[float], data_namings: list[str], white_idx: int, black_idx: int) -> None:
    white_data_average = interpolate(average(all_data[white_idx]))
    black_data_average = interpolate(average(all_data[black_idx]))

    for idx in range(len(all_data)):
        one_file_data = average(all_data[idx])

        if not (idx == white_idx or idx == black_idx):
            one_file_data = normalize(one_file_data, white_data_average, black_data_average)  # ToDO: interpolated?

        one_file_data = median(one_file_data)
        one_file_data = interpolate(one_file_data)

        plot_name = data_namings[idx]
        normalized = True
        if idx == white_idx:
            plot_name = 'White - ' + plot_name
            normalized = False
        elif idx == black_idx:
            plot_name = 'Black - ' + plot_name
            normalized = False
        
        plot(lambda_data, one_file_data, plot_name, normalized)


if __name__ == "__main__":
    # Black and white indexes
    WHITE_IDX = 0
    BLACK_IDX = 2

    # Getting files *.lab in current working directory

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

    analyze(all_data, lambda_data, filenames, WHITE_IDX, BLACK_IDX)

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

# ToDO
# def get_lambda(filename: str) -> list:
#     with open(filename, 'r', encoding='utf16') as f:
#         total_str = ''.join(f.readlines())
#
#     while total_str.find('vecteur') != -1:
#         # Skip string to next vector
#         total_str = total_str[total_str.find('vecteur') + len('vecteur]\n\t'):]
#         total_str = total_str[total_str.find('\n\t') + 2:]
#
#         # Get name
#         name = total_str[:total_str.find('\n')]
#         name = name[name.find('"') + 1:-1]
#
#         # Getting data
#         if name == 'λ':
#             total_str = total_str[total_str.find('points'):]
#             data = total_str[:total_str.find("}") + 1]
#             data = data[data.find("{") + 1:-1].replace("\t", " ").replace("\n", " ").split()
#             start = np.where(np.array(data) >= 400)[0][0]
#             end = np.where(np.array(data) <= 700)[0][-1] + 1
#             return data
#
#     raise NameError('Cannot find lambda in file')


def average(arrays):
    result = []
    for i in range(len(arrays[0])):
        count = len(arrays)
        sum = 0
        for array in arrays:
            sum += array[i]
        av = sum / count  # ToDo: ????? деление?
        result.append(av)
    return result


def median(array):
    result = []
    for i in range(2, len(array) - 2):
        result.append(np.median(array[i-2:i+3]))

    return result


def interpolate(array):
    return savgol_filter(array, 5, 3)


def plot(x, y):
    # ToDo: add x
    matplotlib.pyplot.plot(x, y)
    matplotlib.pyplot.show()


def analyze(all_data, lambda_data):
    for one_file_data in all_data:
        one_file_data = average(one_file_data)
        one_file_data = median(one_file_data)
        one_file_data = interpolate(one_file_data)
        plot(lambda_data, one_file_data)


if __name__ == "__main__":
    # Getting files *.lab in current working directory

    only_files = [f for f in listdir(getcwd()) if isfile(join(getcwd(), f))]
    files = list(filter(lambda x: re.match(r'.*\.lab$', x), only_files))

    # ToDO
    # writer = pandas.ExcelWriter('output.xlsx', engine='xlsxwriter')

    all_data = []
    filenames = []

    naming_pattern = re.compile(r'([0-9]+)\.lab$')
    lambda_data = []

    for i in sorted(files, key=lambda x: naming_pattern.search(x).group().replace('.lab', '')):
        data_array, lambda_data = get_all_vectors_from_file(i)
        all_data.append(data_array)
        filenames.append(naming_pattern.search(i).group().replace('.lab', ''))

    # lambda_data = get_lambda(filenames[0] + '.lab')

    # Cut lambda array due to using median
    lambda_data = lambda_data[2: len(lambda_data) - 2]

    analyze(all_data, lambda_data)

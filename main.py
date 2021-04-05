import pandas
import re
from os import listdir, getcwd
from os.path import isfile, join


# Getting vectors from file
def filter_vectors(filename: str, naming_pattern=r"^\[[0-9]*\]-.*-[0-9]*_E$") -> dict:
    """
    Filtering data from *.lab
    :param filename: name of the file
    :param naming_pattern: which vectors not to ignore
    :return: Array: [ [ 'Name', Array: [ 'Data' ] ] ]
    """

    with open(filename, 'r', encoding='utf16') as f:
        total_str = ''.join(f.readlines())

    lambda_was_already = False
    dictionary = dict()

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

        dictionary[name] = list(map(float, data))

    return dictionary


def export_data_from_lab_to_xlsx(pandas_excel_writer: pandas.ExcelWriter, dictionaries_of_values: list[dict], sheet_names: list[str]) -> None:
    k = 0
    for dictionary in dictionaries_of_values:
        name = sheet_names[k]
        k += 1

        # Write lambda first
        pandas_data_frame = pandas.DataFrame(dictionary)

        # Reorder columns
        lambda_index = pandas_data_frame.columns.tolist().index('λ')
        columns = [pandas_data_frame.columns.tolist()[lambda_index]] + pandas_data_frame.columns.tolist()[:lambda_index] + pandas_data_frame.columns.tolist()[lambda_index+1:]

        pandas_data_frame = pandas_data_frame.reindex(columns, axis=1)

        # Filtering lambda < 400 and lambda > 700
        pandas_data_frame = pandas_data_frame[400 <= pandas_data_frame['λ']]
        pandas_data_frame = pandas_data_frame[700 >= pandas_data_frame['λ']]

        pandas_data_frame.to_excel(excel_writer=pandas_excel_writer, sheet_name=name, index=False)

    pandas_excel_writer.save()


if __name__ == "__main__":
    # Getting files *.lab in current working directory

    only_files = [f for f in listdir(getcwd()) if isfile(join(getcwd(), f))]
    files = list(filter(lambda x: re.match(r'.*\.lab$', x), only_files))

    writer = pandas.ExcelWriter('output.xlsx', engine='xlsxwriter')

    all_data = []
    filenames = []

    naming_pattern = re.compile(r'([0-9]+)\.lab$')

    for i in sorted(files, key=lambda x: naming_pattern.search(x).group().replace('.lab', '')):
        all_data.append(filter_vectors(i))
        filenames.append(naming_pattern.search(i).group().replace('.lab', ''))

    export_data_from_lab_to_xlsx(writer, all_data, filenames)

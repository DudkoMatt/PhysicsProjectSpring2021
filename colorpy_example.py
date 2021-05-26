import colorpy.ciexyz
import colorpy.colormodels

from main import *
# del analyze


# def analyze(all_data: list[list[list[float]]], lambda_data: list[float], data_namings: list[str], white_idx: int,
#             black_idx: int) -> list[list[float]]:
#     """
#     Main function for analyzing, filtering, interpolating, plotting data
#     :param all_data: several samples data 3-d array (different samples -> several measurements -> measurement)
#     :param lambda_data: lambda values for which measurements are taken
#     :param data_namings: names for the samples
#     :param white_idx: white data index in all_data list
#     :param black_idx: black data index in all_data list
#     :return: data after calculations
#     """
#
#     white_data_average = average(all_data[white_idx])
#     black_data_average = average(all_data[black_idx])
#
#     result = []
#
#     for idx in range(len(all_data)):
#         # one_file_data = interpolate(all_data[idx])
#         one_file_data = average(all_data[idx])
#         # one_file_data = average(one_file_data)
#
#         one_file_data = median(one_file_data)
#
#         # if not (idx == white_idx or idx == black_idx):
#         #     one_file_data = normalize(one_file_data, white_data_average, black_data_average)  # ToDO: interpolated?
#
#         # one_file_data = interpolate(one_file_data)
#
#         plot_name = data_namings[idx]
#         normalized = True
#         if idx == white_idx:
#             plot_name = 'White - ' + plot_name
#             normalized = False
#         elif idx == black_idx:
#             plot_name = 'Black - ' + plot_name
#             normalized = False
#
#         result.append(one_file_data)
#         # plot(lambda_data, one_file_data, plot_name, normalized)
#
#     return result

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
        # xyz_coord = colorpy.ciexyz.xyz_from_spectrum(np.array([[lambda_data[i], vector[i]] for i in range(len(vector))]))
        print("{:.5f} {:.5f} {:.5f}".format(*xyz_coord), filenames[k], file=xyz_file)
        k += 1

        calculated_colors.append(colorpy.colormodels.rgb_from_xyz(xyz_coord))
    if platform.system() == "Windows":
        output_file = open('output.txt', 'w')
        for calculated_color in calculated_colors:
            print("{} {} {}".format(*list(map(lambda _x: max(min(round(_x * 255), 255), 0), calculated_color))),
                  file=output_file)

        output_file.flush()
        output_file.close()
    else:
        print(*calculated_colors, sep='\n',
              file=open('/Users/igorklyuzev/ITMO/4_семестр/Физика/PhysicsProjectSpring2021/output.txt',
                        'w'))  # ToDO: wrong answers -> negative rgb coordinates

    k = 0
    for color in calculated_colors:
        img = Image.new('RGB', (1000, 1000), (
        max(min(round(color[0] * 255), 255), 0), max(min(round(color[1] * 255), 255), 0),
        max(min(round(color[2] * 255), 255), 0)))
        img.save("./images/colorpy/{}.png".format(filenames[k]))
        k += 1

    xyz_file.flush()
    xyz_file.close()


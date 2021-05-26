from main import *

def mix_it(first_color: list[float], second_color: list[float], first_color_part, second_color_part) -> list[float]:
    all = first_color_part + second_color_part
    first = first_color_part / all
    second = second_color_part / all
    result = [.0] * len(first_color)
    assert len(first_color) == len(second_color)
    for i in range(0, len(first_color)):
        mixed = (first_color[i] ** first) * (second_color ** second)
        result.append(mixed)
    return result

# ToDO:
'''
1. передать в mix_it интерполированные(или нет) данные со спектрометра
2. получить смешанный результат
3. перевести xyz, потом в rgb и вывести смешанные цвет
4. непонятно как подбирать это по спектру
'''
import math

# this function originates from my AMULATOR_offline
# https://github.com/LangeTo/AMULATOR_offline/blob/main/couplex_functions.py


def add_pos_par(well, doubles, df):
    doubles = doubles + df[df["Well"] == well]["Count categories"].tolist()[0]

    return doubles


# https://www.knowledgehut.com/blog/programming/python-rounding-numbers
def round_up(n, decimals=0):
    """
    Round up a number to a specified number of decimal places.

    Parameters:
    n (float): The number to be rounded up.
    decimals (int): The number of decimal places to round up to. Default is 0.

    Returns:
    float: The number rounded up to the specified number of decimal places.

    Example:
    >>> round_up(2.123, 2)
    2.13
    >>> round_up(2.125, 2)
    2.13
    >>> round_up(2.125, 0)
    3.0
    """

    multiplier = 10**decimals
    return math.ceil(n * multiplier) / multiplier

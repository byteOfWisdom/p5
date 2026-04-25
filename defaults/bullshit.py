import locale
import matplotlib.pyplot as plt


def this_is_fucking_stupid_no_one_actually_gives_a_fuck():
    # Set to German locale to get comma decimal separater
    locale.setlocale(locale.LC_NUMERIC, 'de_DE.UTF-8')

    plt.rcdefaults()

    # Tell matplotlib to use the locale we set above
    plt.rcParams['axes.formatter.use_locale'] = True

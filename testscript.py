import matplotlib.pyplot as plt

def plots(row_number, col_number, width, height, *args, suptitle=None,**kwargs):
    if len(args) != row_number * col_number:
        raise ValueError("Le nombre d'arguments fournis ne correspond pas au nombre de sous-graphiques attendus.")
    
    fig, axs = plt.subplots(row_number, col_number, figsize=(width, height))

    for i, plot_call in enumerate(args):
        col = i // row_number
        row = i % row_number
        plot_call(axs[row, col])

    for ax in axs.flat:
        ax.set(**kwargs)

    plt.suptitle(suptitle)

    plt.tight_layout()
    plt.show()
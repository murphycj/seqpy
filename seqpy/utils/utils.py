import pandas

def sort_by_row_variance(data,descending=True):
    """
    Sort the matrix by variance along the rows
    """

    index = data.index

    variance = data.var(axis=1)

    sort_index = sorted(range(len(variance)), key=lambda k: variance[k], reverse=descending)

    index = index[sort_index]
    data = data.ix[index]

    return data

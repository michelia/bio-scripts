def two_dimensional_array_transformation(arr):
    """
    In [1]: arr = [[1, 2],[3, 4],[5, 6]]

    In [2]: [[r[col] for r in arr] for col in range(len(arr[0]))]
    Out[2]: [[1, 3, 5], [2, 4, 6]]

    In [3]: zip(*arr)
    Out[3]: [(1, 3, 5), (2, 4, 6)]

    In [4]: map(list, zip(*arr))
    Out[4]: [[1, 3, 5], [2, 4, 6]]
    ---------------------------------------------------------------
    In [1]: arr = [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]]

    In [2]: [[r[col] for r in arr] for col in range(len(arr[0]))]
    Out[2]: [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]

    In [3]: zip(*arr)
    Out[3]: [(1, 4, 7, 10), (2, 5, 8, 11), (3, 6, 9, 12)]

    In [4]: map(list, zip(*arr))
    Out[4]: [[1, 4, 7, 10], [2, 5, 8, 11], [3, 6, 9, 12]]
    """

    # return a list of sub_list, but the lower speed than below
    # return [[r[col] for r in arr] for col in range(len(arr[0]))]

    # return a list of sub_tuple
    return zip(*arr)

    # return a list of sub_list
    # return map(list, zip(*arr))

    # return a iter, same as map(list, zip(*arr))
    # import itertools
    # return map(list, itertools.izip(*arr))

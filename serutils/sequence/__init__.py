
def nicer(res, sep=' '):
    """
    writes resolution number for human beings.
    """
    if not res % 1000000000:
        return str(res)[:-9] + sep + 'Gb'
    if not res % 1000000:
        return str(res)[:-6] + sep + 'Mb'
    if not res % 1000:
        return str(res)[:-3] + sep + 'kb'
    if res == 1:
        return 'bin'
    return str(res) + sep + 'b'

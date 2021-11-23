import gzip

def get_file_handle(file_name, mode):

    if file_name.endswith(".gz"):
        fh = gzip.open(file_name, mode)
    else:
        fh = open(file_name, mode)

    return(fh)
import gzip

def get_file_handle(file_name, mode):

    if file_name.endswith(".gz"):
        handle = gzip.open(file_name, mode)
    else:
        handle = open(file_name, mode)

    return(handle)
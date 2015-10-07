import linecache
import sys

def PrintException():
    '''
    modified from http://stackoverflow.com/questions/14519177/python-exception-handling-line-number
    '''
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print 'Exception raised from {}:{}, "{}":\n{}: {}'.format(filename, lineno, line.strip(), exc_type.__name__, exc_obj)
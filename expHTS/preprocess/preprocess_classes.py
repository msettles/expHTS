# preprocess classes
import os
from itertools import izip_longest


def which(program):
    '''
    mimics 'which' on the command line
    borrowed from
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    '''
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def version_check(appVersion, requiredVersion):
    '''
    Compare version numbers
    assume version number is '.' delimeted
    '''
    for i, j in izip_longest(appVersion.split('.'), requiredVersion.split('.'), fillvalue=None):
        if i is None or j is None:
            return False
        if i < j:
            return False
        if i > j:
            return True
    return True


##############################################################################
# Wrappers for preprocessing apps
##############################################################################
######
# Each app has 4 functions
# app_check - version checker
# app_call - builds the command line string
# app_output - parses the output
# app_run - pulls everything together and runs the app
######
class wrapper_base(object):
    '''
    base class for preprocess application wrappers
    '''
    def app_check(self):
        '''
        verify the application existance and version
        '''
        pass

    def app_call(self):
        '''
        build the call to the app
        '''
        pass

    def app_output(self):
        '''
        parse the output
        '''
        pass

    def app_run(self, args):
        '''
        put it all together
        '''
        pass

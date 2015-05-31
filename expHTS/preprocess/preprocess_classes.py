# preprocess classes
import os
import re
from itertools import izip_longest
from subprocess import Popen
from subprocess import PIPE


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
    def __init__(self):
        pass

    def app_check(self):
        '''
        verify the application existance and version
        '''
        return True

    def app_call(self):
        '''
        build the call to the app
        '''
        return None

    def app_output(self):
        '''
        parse the output
        '''
        return None

    def app_run(self, args):
        '''
        put it all together
        '''
        return 0


class read_processing(wrapper_base):
    '''
    Process the final set of reads at the end of preproccess
    remove final too-short reads, change reads names (SRA, Illumina new to old)
    Count reads, length and quality in PE,SE
    Produce final set of sequence files
    '''
    def __init__(self, inputs, outputs, obj, args):
        self.inputs = inputs
        self.output = outputs
        self.object = obj
        self.args = args


class sickle_processing(wrapper_base):
    '''
    Run sickle, to trim for Quality and polyA/T
    '''
    def __init__(self, inputs, outputs, obj, args):
        self.inputs = inputs
        self.output = outputs
        self.object = obj
        self.path = self.app_check("sickle", "1.33")
        self.args = {'-l': 150, '-q': 24, '-t': 'illumina'}
        for arg in args:
            self.args[arg] = arg

    def app_check(self, name, version):
        sickle_path = which(name)
        if sickle_path is None:
            raise Exception('sickle not found on the path')
        try:
            res = Popen([sickle_path, '--version'], stdout=PIPE)
        except OSError as e:
            print "OS error({0}): {1}".format(e.errno, e.strerror)
            raise
        app_version = res.stdout.readline().rstrip()
        app_version = re.split(' +', app_version)[2]
        if version_check(app_version, version) is False:
            raise Exception('require sickle version %s (system has version %s' % (version, app_version))
        return sickle_path

    def app_call(self):
        return [self.app_pe_call(), self.app_se_call()]

    def app_pe_call(self):
        try:
            sickle_pe_call = [self.path, 'pe']
            for key, value in self.args.iteritems():
                temp = [key, value]
                sickle_pe_call.extend(temp)
            sickle_pe_call.extend(
                ['-c', self.inputs['interleaved'],
                 '-M', self.output['stdout']])
            return sickle_pe_call
        except KeyError:
            return None
        except Exception as e:
            print e.message, e.args
            raise

    def app_se_call(self):
        try:
            sickle_se_call = [self.path, 'se']
            for key, value in self.args.iteritems():
                temp = [key, value]
                sickle_se_call.extend(temp)
            sickle_se_call.extend(
                ['-f', self.inputs['se'],
                 '-o', self.output['stdout']])
            return sickle_se_call
        except KeyError:
            return None
        except Exception as e:
            print e.message, e.args
            raise

    def app_output(self, output_lines):
        sickle_output = {"Reads": 0,
                         "PE1_discarded": 0,
                         "PE1_left-trim": 0,
                         "PE1_right-trim": 0,
                         "PE2_discarded": 0,
                         "PE2_left-trim": 0,
                         "PE2_right-trim": 0,
                         "Both_discarded": 0,
                         "Pairs_kept": 0,
                         "Single_kept": 0,
                         "Pairs_kept_percent": 0,
                         "Reads_kept_percent": 0,
                         "PE1_Poly-AT": 0,
                         "PE2_Poly-AT": 0
                         }
        for line in output_lines:
            line = line.split(": ")
            if (line[0] == "Total input FastQ records" or line[0] == "Total FastQ records"):
                line = re.sub(',|\(|\)', '', line)
                sickle_output['Reads'] = line.split(" ")[0]
            elif line[0] == "FastQ paired records discarded":
                pass
            elif line[0] == "FastQ single records discarded" or line[0] == "FastQ records discarded":
                pass
            elif line[0] == "Poly AT tail PE1" or line[0] == "Poly AT tails":
                pass
            elif line[0] == "Poly AT tail PE2":
                pass
            elif line[0] == "PE1 Base pairs left removed" or line[0] == "Base pairs left removed":
                pass
            elif line[0] == "PE1 Base pairs right removed" or line[0] == "Base pairs right removed":
                pass
            elif line[0] == "PE2 Base pairs left removed":
                pass
            elif line[0] == "PE2 Base pairs right removed":
                pass
            else:
                pass
        return sickle_output



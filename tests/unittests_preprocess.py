import unittest
import expHTS


class preprocess_classes_tests(unittest.TestCase):
    def test_version_check_true(self):
        self.assertTrue(expHTS.preprocess.preprocess_classes.version_check("1.2.0", "1.1.0"))
        self.assertTrue(expHTS.preprocess.preprocess_classes.version_check("1.2.0", "1.1"))

    def test_version_check_false(self):
        self.assertFalse(expHTS.preprocess.preprocess_classes.version_check("1.1.0", "1.2.0"))
        self.assertFalse(expHTS.preprocess.preprocess_classes.version_check("1.1.0", "1.1"))


class sampleSheet_classes_tests(unittest.TestCase):
    def test_reading_sample_sheet(self):
        ss = expHTS.samplesheet.sampleSheet('samples.txt')
        self.assertTrue(ss.getSampleCount() is 4)       


if __name__ == '__main__':
    unittest.main()

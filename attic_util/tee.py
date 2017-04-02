import sys

class Tee(object):
    def __init__(self, fptr):
        self.file = fptr

    def __enter__(self):
        self.stdout = sys.stdout
        sys.stdout = self

    def __exit__(self, exception_type, exception_value, traceback):
        sys.stdout = self.stdout

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)

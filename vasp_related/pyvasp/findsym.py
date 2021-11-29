class findsym:
    def __init__(self):
        self._file_ = 'findsym.log'
        self.lines = ''

    def read(self, file=None):
        if file is not None:
            try:
                file.pos
            except AttributeError:
                self._file_ = open(file)
            else:
                self._file_ = file
                self._file_.seek(0)
        else:
            self._file_ = open(self._file_)

        self.lines = self._file_.readlines()

        return self
    
    def get_spacegroup(self):
        for line in self.lines:
            if line.find('Space Group') > -1:
                line_split = line.split()
                break
        return (line_split[2:])


    def get_latttice_parameters(self):
        for i, line in enumerate(self.lines):
            if line.find('Values of') > -1:
                latpar = self.lines[i+1].split()
                break
        return latpar
            

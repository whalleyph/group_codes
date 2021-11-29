import re

class EigenvalReader:
    def __init__(self):
        self.listeners = {}

    def add_listener(self, event_name, callback):
        if (not self.listeners.has_key(event_name)):
            self.listeners[event_name] = []
        self.listeners[event_name].append(callback)

    def __raise__(self, event_name, args):
        if not self.listeners.has_key(event_name):
            return

        callbacks = self.listeners[event_name]
        if callbacks != None:
            for c in callbacks:
                c(self, args)

    def read_file(self, eigen_file):
        f = open(eigen_file)
        for i in range(0,5):
            f.readline()

        throwaway, kpoint_count, band_count = f.readline().split()
        self.__raise__('ON_KPOINT_COUNT', kpoint_count)
        self.__raise__('ON_BAND_COUNT', band_count)

        # Blank Line
        f.readline()

        # End of File (EOF)
        eof = False
        while (not eof):
            line = f.readline()
            if line == None:
                eof = True
                continue
            try:
                x,y,z,weight = line.split()
                self.__raise__('ON_KPOINT', { "x": x, "y": y, "z": z, "weight": weight })
            except:
                print('failed to split line: ' + line)

            for i in range(0, int(band_count)):
                line = f.readline()
                try:
                    id, energy = line.split()
                    self.__raise__('ON_BAND', { 'id': id, 'energy': float(energy) })
                except:
                    print('failed to get band info from: ' + line)

            line = f.readline()
            if line == '':
                eof = True

class EigenvalParser:
    def parse(self, filename):
        self.kpoints = []
        reader = EigenvalReader()
        reader.add_listener('ON_KPOINT', self.__on_kpoint__)
        reader.add_listener('ON_BAND', lambda reader, band: self.kpoints[-1]['bands'].append(band))

        reader.read_file(filename)

        return self.kpoints

    def __on_kpoint__(self, reader, kpoint):
        kpoint['bands'] = []
        self.kpoints.append(kpoint)

def get_bands(kpoints):
    bands = [[] for b in range(0,len(kpoints[0]['bands']))]

    for kp in kpoints:
        kp_bands = kp['bands']
        for kp_band in kp_bands:
            bands[int(kp_band['id']) - 1].append(kp_band['energy'])

    return bands

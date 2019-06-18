from typing import List, NamedTuple
import gfapy
import pickle
import subprocess
import io

class GFA:
    def __init__(gfa: gfapy.Gfa):
        self.gfa = gfa
    
    def load_from_pickle(file: str):
        self.gfa = pickle.load(file)

    def load_form_xg(file: str, xg_bin: str):
        #self.gfa = XGWrapper.load(file, xg_bin)
        gfa = subprocess.check_output([xg_bin, "-i", file, "--gfa-out"])
        process = subprocess.Popen([xg_bin, "-i", file, "--gfa-out"], stdout=subprocess.PIPE)
        with io.open(process.stdout.fileno(), closefd=False) as stream:
            [self.gfa.add_line(line) for line in stream]
        process.wait()
        if process.returncode != 0:
            raise OSError()

    def load_from_gfa(file: str):
        self.gfa = gfapy.Gfa.from_file(file) 

    def save_as_pickle(self):
        pickle.dump(self.gfa, file)

    def save_as_xg(self, file: str):
        process = subprocess.Popen([xg_bin, "-o", file, "-g", "-"], stdin=subprocess.PIPE)
        output = process.communicate(self.gfa.to_gfa1_s())
        return output

    def save_as_gfa(self, file:str):
        self.gfa.to_file(file)

'''
class XGWrapper:
    @staticmethod
    def save(gfa):
        pass
    
    @staticmethod
    def load(gfa):
        pass

class GraphStack:
    def __init__(graphs: List[Graph]):
        self.graphs = graphs
'''

if __name__ == "__main__":
    location_of_xg = sys.argv[0]



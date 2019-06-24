from typing import List, NamedTuple
import gfapy
import pickle
import subprocess
import io
import os
import tempfile

class GFA:
    def __init__(self, gfa: gfapy.Gfa):
        self.gfa = gfa
    
#    @classmethod
#    def load_from_pickle(cls, file: str):
#        graph = pickle.load(open(file, 'rb'))
#        return graph

    @classmethod
    def load_form_xg(cls, file: str, xg_bin: str):
        gfa = gfapy.Gfa()
        process = subprocess.Popen([xg_bin, "-i", file, "--gfa-out"], stdout=subprocess.PIPE)
        with io.open(process.stdout.fileno(), closefd=False) as stream:
            [gfa.add_line(line) for line in stream]
        process.wait()
        if process.returncode != 0:
            raise OSError()
        graph = cls(gfa)
        process.stdout.close()
        return graph

    @classmethod
    def load_from_gfa(cls, file: str):
        gfa = gfapy.Gfa.from_file(file) 
        graph = cls(gfa)
        return graph

#    def save_as_pickle(self, outfile: str):
#        with open(outfile, 'wb') as pickle_file:
#            pickle.dump(self.gfa, pickle_file, protocol=2)

    def save_as_xg(self, file: str, xg_bin:str):
        with tempfile.NamedTemporaryFile(mode="w+t", delete=False) as f:
            f.write(self.gfa.to_gfa1_s())

        process = subprocess.check_output([xg_bin, "-o", file, "-g", f.name])
        os.remove(f.name)

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



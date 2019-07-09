from src.graph import *

import dataclasses

@dataclasses.dataclass
class Profile:
    node: NodeIndex
    paths: List[Path]
    duplicate: int = 0


class DAGify:
    def __init__(self, paths: List[Path], nodes = {}):
        """

        :type paths: List[Path]
        """
        self.paths = paths
        self.nodes = nodes
        self.profile = []

    # def random_search_to_minimize_node_replication(self):


    def recursive_merge(self, primary_path_index: int = 0):
        profile = []
        for node_index in self.paths[primary_path_index].nodes:
            profile.append(Profile(node_index, [self.paths[primary_path_index]], 0))
        for i, path in enumerate(self.paths):
            if i == primary_path_index:
                continue
            profile = self.lcs(profile, path)
        return profile

    def lcs(self, s1: List[Profile], s2: Path) -> List[Profile]:
        n, m = len(s1), len(s2.nodes)
        dp = [[0] * (m+1) for _ in range(n+1)]

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if s1[i-1].node == s2.nodes[j-1]:
                    dp[i][j] = dp[i - 1][j - 1] + 1
                else:
                    dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
        i, j = n, m
        index = []
        prev = set()

        while i > 0 and j > 0:
            if s1[i-1].node == s2.nodes[j-1]:
                prev_paths = s1[i-1].paths
                prev_paths.append(s2)
                index.append(Profile(s1[i-1].node, prev_paths, s1[i-1].node.node.index in prev))
                prev.add(s1[i-1].node.node.index)
                i -= 1
                j -= 1
            elif dp[i-1][j] > dp[i][j-1]:
                prev_paths = s1[i-1].paths
                index.append(Profile(s1[i-1].node, prev_paths, s1[i-1].node.node.index in prev))
                prev.add(s1[i-1].node.node.index)
                i -= 1
            else:
                index.append(Profile(s2.nodes[j-1], [s2], False))
                prev.add(s2.nodes[j-1].node.index)
                j -= 1

        while i > 0:
            prev_paths = s1[i - 1].paths
            index.append(Profile(s1[i - 1].node, prev_paths, s1[i - 1].node.node.index in prev))
            prev.add(s1[i - 1].node.node.index)
            i -= 1

        while j > 0:
            prev.add(s2.nodes[j - 1])
            index.append(Profile(s2.nodes[j - 1], [s2], False))
            j -= 1

        index.reverse()
        self.profile = index

        return index

    def to_graph(self):
        factory_input = []
        current_slice = Slice([])
        print(self.profile)
        for prof in self.profile:
            paths = [x.name for x in prof.paths]
            if len(prof.paths) == len(self.paths):
                if len(current_slice.nodes) > 0:
                    factory_input.append(current_slice)
                factory_input.append(Slice([Node(prof.node.node.seq, paths, prof.node.node.index)]))
                current_slice = Slice([])
            else:
                all_set = set()
                for items in [x.paths for x in current_slice.nodes]:
                    all_set = all_set | items
                if set(prof.paths) & all_set != set():
                    if len(current_slice.nodes) > 0:
                        current_slice.add_node(Node("", set([x.name for x in self.paths]) - all_set))
                        factory_input.append(current_slice)
                    current_slice = Slice([Node(prof.node.node.seq, paths, prof.node.node.index)])
                else:
                    current_slice.add_node(Node(prof.node.node.seq, paths, prof.node.node.index))

        base_graph = Graph.load_from_slices(factory_input)
        print(factory_input)
        return base_graph

    def merge(A: List[NodeIndex], B: List[NodeIndex]):
        pos, merged = [], []
        pi, pj, prev = 0, 0, set()
        for i in range(len(A)):
            for j in range(len(B)):
                if pi <= i and pj <= j and A[i] == B[j]:
                    curr = set()
                    while pi < i:
                        curr.add(A[pi])
                        pos.append( (pi, -1, A[pi] in prev) )
                        merged.append(A[pi])
                        pi += 1
                    while pj < j:
                        curr.add(B[pj])
                        pos.append( (-1, pj, B[pj] in prev) )
                        merged.append(B[pj])
                        pj += 1
                    if i == pi and j == pj:
                        pos.append((i, j, False))
                        merged.append(A[i])
                        pi += 1
                        pj += 1
                    prev |= curr
        while pi < len(A):
            pos.append( (pi, -1, A[pi] in prev) )
            merged.append(A[pi])
            pi += 1
        while pj < len(B):
            pos.append( (-1, pj, B[pj] in prev) )
            merged.append(B[pj])
            pj += 1
        return pos, merged
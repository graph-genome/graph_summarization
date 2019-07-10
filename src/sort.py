from src.graph import *

import dataclasses

@dataclasses.dataclass
class Profile:
    node: NodeIndex
    paths: List[Path]
    candidate_paths: set()
    duplicate: bool = False


class DAGify:
    def __init__(self, paths: List[Path], nodes = {}):
        """

        :type paths: List[Path]
        """
        self.paths = paths
        self.nodes = nodes
        self.profile = []

    def search_for_minimizing_replications(self) -> (List[Profile], int):
        min_rep = len(self.nodes)
        profile = []
        for i, _ in enumerate(self.paths):
            profile_candidate = self.recursive_merge(i)
            if min_rep > sum([x for x in profile if x]):
                min_rep = sum([x for x in profile if x])
                profile = profile_candidate
        return profile, min_rep

    def recursive_merge(self, primary_path_index: int = 0) -> List[Profile]:
        profile = []
        for node_index in self.paths[primary_path_index].nodes:
            profile.append(Profile(node_index, [self.paths[primary_path_index]], {self.paths[primary_path_index].name}, 0))
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
        candidate_path_flag = False

        while i > 0 and j > 0:
            if s1[i-1].node == s2.nodes[j-1]:
                prev_paths = s1[i-1].paths
                prev_paths.append(s2)
                candidate_paths = s1[i-1].candidate_paths
                candidate_paths.add(s2.name)
                candidate_path_flag = True

                index.append(Profile(s1[i-1].node, prev_paths, candidate_paths, s1[i-1].node.node.index in prev))
                prev.add(s1[i-1].node.node.index)
                i -= 1
                j -= 1
            elif dp[i-1][j] > dp[i][j-1]:
                prev_paths = s1[i-1].paths
                candidate_paths = s1[i-1].candidate_paths
                if candidate_path_flag:
                    candidate_paths.add(s2.name)
                index.append(Profile(s1[i-1].node, prev_paths, candidate_paths, s1[i-1].node.node.index in prev))
                prev.add(s1[i-1].node.node.index)
                i -= 1
            else:
                candidate_paths = {s2.name}
                if s1[i]:
                    candidate_paths |= s1[i].candidate_paths
                if s1[i-1]:
                    candidate_paths |= s1[i-1].candidate_paths
                index.append(Profile(s2.nodes[j-1], [s2], candidate_paths, False))
                prev.add(s2.nodes[j-1].node.index)
                j -= 1

        while i > 0:
            prev_paths = s1[i - 1].paths
            prev_candidates = s1[i-1].candidate_paths
            index.append(Profile(s1[i - 1].node, prev_paths, prev_candidates, s1[i - 1].node.node.index in prev))
            prev.add(s1[i - 1].node.node.index)
            i -= 1

        while j > 0:
            prev.add(s2.nodes[j - 1])
            index.append(Profile(s2.nodes[j - 1], [s2], {s2.name}, False))
            j -= 1

        index.reverse()
        self.profile = index

        return index

    def to_graph(self):
        factory_input = []
        current_slice = Slice([])
        # print(self.profile)
        for prof in self.profile:
            paths = [x.name for x in prof.paths]
            if len(prof.paths) == len(prof.candidate_paths):
                if len(current_slice.nodes) > 0:
                    factory_input.append(current_slice)
                factory_input.append(Slice([Node(prof.node.node.seq, paths, prof.node.node.index)]))
                current_slice = Slice([])
            else:
                all_set = set()
                for items in [x.paths for x in current_slice.nodes]:
                    all_set |= items
                # print(all_set, prof.candidate_paths, prof.paths, set([x.name for x in prof.paths]) & all_set)
                if set([x.name for x in prof.paths])  & all_set != set():
                    if len(current_slice.nodes) > 0:
                        if prof.candidate_paths - all_set != set():
                            current_slice.add_node(Node("", prof.candidate_paths - all_set))
                        factory_input.append(current_slice)
                    current_slice = Slice([Node(prof.node.node.seq, paths, prof.node.node.index)])
                else:
                    current_slice.add_node(Node(prof.node.node.seq, paths, prof.node.node.index))

        base_graph = Graph.load_from_slices(factory_input)
        # print(factory_input)
        return base_graph

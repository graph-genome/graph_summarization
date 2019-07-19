from src.graph import *

import dataclasses

@dataclasses.dataclass
class Profile:
    node: NodeTraversal
    forward_paths: List[Path]
    backward_paths: List[Path]
    candidate_paths: set()
    duplicate: bool = False

    def __repr__(self):
        return "[" + str(self.node.node) + str(self.forward_paths) + str(self.backward_paths) + ":" + str(self.candidate_paths) + "]"

class DAGify:
    def __init__(self, paths: List[Path], nodes={}):
        """

        :type paths: List[Path]
        """
        self.paths = paths
        self.nodes = nodes

    def search_for_minimizing_replications(self) -> (List[Profile], int):
        min_rep = sys.maxsize
        profile = []
        for i, _ in enumerate(self.paths):
            profile_candidate = self.recursive_merge(i)
            if min_rep > len([x.duplicate for x in profile_candidate if x.duplicate]):
                min_rep = len([x.duplicate for x in profile_candidate if x.duplicate])
                profile = profile_candidate
        return profile, min_rep

    def recursive_merge(self, primary_path_index: int = 0) -> List[Profile]:
        profile = []
        for node_index in self.paths[primary_path_index].nodes:
            if node_index.strand == "+":
                profile.append(Profile(node_index, [self.paths[primary_path_index]], [], {self.paths[primary_path_index]}, False))
            else:
                profile.append(
                    Profile(node_index, [], [self.paths[primary_path_index]], {self.paths[primary_path_index]}, False))
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
            if s1[i-1].node.node.id == s2.nodes[j-1].node.id:
                prev_fwd_paths = s1[i-1].forward_paths
                prev_bwd_paths = s1[i-1].backward_paths
                if s2.nodes[j-1].strand == "+":
                    prev_fwd_paths.append(s2)
                else:
                    prev_bwd_paths.append(s2)
                candidate_paths = s1[i-1].candidate_paths
                candidate_paths.add(s2)
                candidate_path_flag = True

                index.append(Profile(s1[i-1].node, prev_fwd_paths, prev_bwd_paths, candidate_paths, s1[i-1].node.node.id in prev))
                prev.add(s1[i-1].node.node.id)
                i -= 1
                j -= 1
            elif dp[i-1][j] > dp[i][j-1]:
                prev_fwd_paths = s1[i-1].forward_paths
                prev_bwd_paths = s1[i-1].backward_paths
                candidate_paths = s1[i-1].candidate_paths
                if candidate_path_flag:
                    candidate_paths.add(s2)
                index.append(Profile(s1[i-1].node, prev_fwd_paths, prev_bwd_paths, candidate_paths, s1[i-1].node.node.id in prev))
                prev.add(s1[i-1].node.node.id)
                i -= 1
            else:
                candidate_paths = {s2}
                if s2.nodes[j-1].strand == "+":
                    fwd_paths = [s2]
                    bwd_paths = []
                else:
                    fwd_paths = []
                    bwd_paths = [s2]
                if i > n and s1[i]:
                    candidate_paths |= s1[i].candidate_paths
                if s1[i-1]:
                    candidate_paths |= s1[i-1].candidate_paths
                index.append(Profile(s2.nodes[j-1], fwd_paths, bwd_paths, candidate_paths, s2.nodes[j-1].node.id in prev))
                prev.add(s2.nodes[j-1].node.id)
                j -= 1

        while i > 0:
            prev_fwd_paths = s1[i - 1].forward_paths
            prev_bwd_paths = s1[i - 1].backward_paths
            prev_candidates = s1[i-1].candidate_paths
            index.append(Profile(s1[i - 1].node, prev_fwd_paths, prev_bwd_paths, prev_candidates, s1[i - 1].node.node.id in prev))
            prev.add(s1[i - 1].node.node.id)
            i -= 1

        while j > 0:
            if s2.nodes[j - 1].strand == "+":
                fwd_paths = [s2]
                bwd_paths = []
            else:
                fwd_paths = []
                bwd_paths = [s2]
            prev.add(s2.nodes[j - 1].node.id)
            index.append(Profile(s2.nodes[j - 1], fwd_paths, bwd_paths, {s2}, False))
            j -= 1

        index.reverse()
        # print(index)

        return index

    def to_slices(self, profile: List[Profile]) -> List[Path]:
        factory_input = []
        current_slice = Slice([])
        current_paths = []

        for index, prof in enumerate(profile):
            fwd_paths = prof.forward_paths
            bwd_paths = prof.backward_paths
            all_path_set = set([x for x in current_paths])
            # print(prof, current_slice, current_paths)
            candidate_paths_set = prof.candidate_paths
            if index + 1 != len(profile):
                candidate_paths_set |= profile[index+1].candidate_paths

            if len(fwd_paths) + len(bwd_paths) == len(candidate_paths_set):
                if len(current_slice.nodes) > 0:
                    if prof.candidate_paths - all_path_set != set():
                        current_slice.add_node(NodeTraversal(Node("", prof.candidate_paths - all_path_set)))
                    factory_input.append(current_slice)
                current_slice = Slice([])
                if fwd_paths != []:
                    current_slice.add_node(NodeTraversal(Node(prof.node.node.seq, fwd_paths, prof.node.node.id), "+"))
                if bwd_paths != []:
                    current_slice.add_node(NodeTraversal(Node(prof.node.node.seq, bwd_paths, prof.node.node.id), "-"))
                factory_input.append(current_slice)
                current_slice = Slice([])
                current_paths = []
            else:
                if (set(fwd_paths) | set(bwd_paths) ) & all_path_set != set():
                    if len(current_slice.nodes) > 0:
                        if prof.candidate_paths - all_path_set != set():
                            current_slice.add_node(NodeTraversal(Node("", prof.candidate_paths - all_path_set)))
                        factory_input.append(current_slice)
                    current_slice = Slice([])
                    if fwd_paths != []:
                        current_slice.add_node(NodeTraversal(Node(prof.node.node.seq, fwd_paths, prof.node.node.id), "+"))
                    if bwd_paths != []:
                        current_slice.add_node(NodeTraversal(Node(prof.node.node.seq, bwd_paths, prof.node.node.id), "-"))
                    current_paths = fwd_paths
                    current_paths.extend(bwd_paths)
                else:
                    if fwd_paths != []:
                        current_slice.add_node(NodeTraversal(Node(prof.node.node.seq, fwd_paths, prof.node.node.id), "+"))
                    if bwd_paths != []:
                        current_slice.add_node(NodeTraversal(Node(prof.node.node.seq, bwd_paths, prof.node.node.id), "-"))
                    current_paths.extend(bwd_paths)
                    current_paths.extend(fwd_paths)

        if len(current_slice.nodes) > 0:
            all_path_set = set([x for x in current_paths])
            if profile[-1].candidate_paths - all_path_set != set():
                current_slice.add_node(NodeTraversal(Node("", prof.candidate_paths - all_path_set)))
            factory_input.append(current_slice)
        return factory_input

    def to_graph(self, profile: List[Profile]):
        factory_input = self.to_slices(profile)
        base_graph = SlicedGraph.load_from_slices(factory_input, self.paths)
        return base_graph

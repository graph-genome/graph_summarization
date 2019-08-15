import sys
import dataclasses
from typing import List

from src.graph import NodeTraversal, Path, Slice, Node, SlicedGraph


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
    """
    DAGify accepts a set of paths, and
    """
    def __init__(self, paths: List[Path], nodes={}):
        """
        :type paths: List[Path], nodes: Set[Node]
        """
        self.paths = paths
        self.nodes = nodes

    def generate_profiles_with_minimizing_replications(self) -> (List[Profile], int):
        """
        Generate profiles with minimizing the number of node repication by trying to assign each path as a primary path one by one.
        It returns profiles and whose number of replication.
        :return: Profiles and the number of replicated nodes.
        """
        min_rep = sys.maxsize
        profile = []
        for i, _ in enumerate(self.paths):
            profile_candidate = self.generate_profiles(i)
            if min_rep > len([x.duplicate for x in profile_candidate if x.duplicate]):
                min_rep = len([x.duplicate for x in profile_candidate if x.duplicate])
                profile = profile_candidate
        return profile, min_rep

    def generate_profiles(self, primary_path_index: int = 0) -> List[Profile]:
        """
        Generate profiles of paths.
        :param primary_path_index: an index of paths to assign as the primary path.
        :return: a list of profiles
        """
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
        """
        Compute longest common substrings between a profile and a path, and return as a new profile.
        Longest common substrings are calculated by dynamic programming on a 2-dimensional array, i.e. O(n^2).
        LCS distinguish the alignment between node ids of paths as a match, mismatch or gap.
        A match: paths are stored on one profile.
        A mismatch: paths are stored on separated node's profile.
        A gap: paths are stored as candidate_paths on a gapped node's profile.

        :param s1: a list of profiles to be merged a path.
        :param s2: a path to merge into a list of profiles.
        :return: a list of profiles
        """
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

    def to_slices(self, profiles: List[Profile]) -> List[Slice]:
        """
        Convert profiles to a list of slices by merging adjacent profiles that have disjoint set of paths into one slice.
        :param profiles: a list of profiles.
        :return: a list of slices.
        """
        factory_input = []
        current_slice = Slice([])
        current_paths = []

        for index, prof in enumerate(profiles):
            fwd_paths = prof.forward_paths
            bwd_paths = prof.backward_paths
            all_path_set = set([x for x in current_paths])
            # print(prof, current_slice, current_paths)
            candidate_paths_set = prof.candidate_paths
            if index + 1 != len(profiles):
                candidate_paths_set |= profiles[index + 1].candidate_paths

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
            if profiles[-1].candidate_paths - all_path_set != set():
                current_slice.add_node(NodeTraversal(Node("", profile[-1].candidate_paths - all_path_set)))
            factory_input.append(current_slice)
        return factory_input

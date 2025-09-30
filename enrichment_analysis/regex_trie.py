from typing import List, Dict

class RegexTrieNode:
    def __init__(self):
        # Edges: map from bitmask of allowed letters → child node
        self.children: Dict[int, 'RegexTrieNode'] = {}
        self.proteases: List[str] = []

class RegexTrie:
    def __init__(self, alphabet: List[str]):
        self.root = RegexTrieNode()
        self.alphabet = alphabet
        self.char_to_bit = {ch: 1 << i for i, ch in enumerate(alphabet)}
        self.full_mask = (1 << len(alphabet)) - 1

    def _expand_set(self, token_list: List[str]) -> int:
        '''Convert a token into a bitmask of allowed letter.'''
        if token_list == ["X"]:
            return self.full_mask
        elif all(tok.startswith("!") for tok in token_list):
            forbidden = {tok[1:] for tok in token_list}
            mask = self.full_mask
            for f in forbidden:
                mask &= ~self.char_to_bit[f]
            return mask
        else:
            mask = 0
            for ch in token_list:
                mask |= self.char_to_bit[ch]
            return mask

    def insert(self, regex: List[List[str]], protease_name: str):
        if regex != [["X"]]*8:
            node = self.root
            for token_list in regex:
                allowed_mask = self._expand_set(token_list)
                if allowed_mask not in node.children:
                    node.children[allowed_mask] = RegexTrieNode()
                node = node.children[allowed_mask]
            if protease_name not in node.proteases:
                node.proteases.append(protease_name)

    def match(self, word: str) -> List[str]:
        matches = []

        def dfs(node: RegexTrieNode, idx: int):
            nonlocal matches
            if idx == len(word):
                matches.extend(node.proteases)
                return
            ch = word[idx]
            if ch not in self.char_to_bit:
                return
            bit = self.char_to_bit[ch]
            for mask, child in node.children.items():
                if mask & bit:
                    dfs(child, idx + 1)

        dfs(self.root, 0)
        return matches
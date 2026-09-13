
import itertools
import json
from pathlib import Path


wenjian = Path(__file__).resolve().parent / "five_point_cliques.json"
xuyao = (1, *range(4, 20))


def bian(a, b):
    if a < b:
        return a, b
    return b, a


def dushuju():
    d = json.loads(wenjian.read_text(encoding="utf-8"))
    if set(d) != {"near_one_factors", "cliques"}:
        raise ValueError("unexpected fields in five_point_cliques.json")

    fen = []
    for x in d["near_one_factors"]:
        fen.append(frozenset(bian(*y) for y in x))
    fen = tuple(fen)

    tuan = {}
    for k, v in d["cliques"].items():
        tuan[int(k)] = tuple(tuple(x) for x in v)
    if tuple(sorted(tuan)) != xuyao:
        raise ValueError("the required clique orders are 1 and 4 through 19")
    return fen, tuan


def chafen(fen):
    guding = (
        frozenset({(1, 4), (2, 3)}),
        frozenset({(0, 2), (3, 4)}),
        frozenset({(0, 4), (1, 3)}),
        frozenset({(0, 1), (2, 4)}),
        frozenset({(0, 3), (1, 2)}),
    )
    if fen != guding:
        raise ValueError("the near-one-factors differ from those in the paper")

    quan = set()
    for i in range(5):
        for j in range(i + 1, 5):
            quan.add((i, j))

    yongguo = []
    for i in range(5):
        f = fen[i]
        dian = set()
        for x in f:
            dian.update(x)
            yongguo.append(x)
        if dian != set(range(5)) - {i}:
            raise ValueError(f"factor {i} is not a matching on the other vertices")
    if len(yongguo) != 10 or set(yongguo) != quan:
        raise ValueError("the near-one-factors do not partition E(K_5)")
    if len(set(yongguo)) != len(yongguo):
        raise ValueError("an edge occurs in two near-one-factors")


def jiantu(fen):
    """Join two permutations when a coordinate pair lies in its factor."""
    pai = list(itertools.permutations(range(5)))
    lin = [set() for x in pai]
    bianshu = 0
    for j in range(len(pai)):
        for i in range(j):
            lian = False
            for k in range(5):
                if bian(pai[i][k], pai[j][k]) in fen[k]:
                    lian = True
                    break
            if lian:
                lin[i].add(j)
                lin[j].add(i)
                bianshu += 1
    return pai, lin, bianshu


def chatuan(jie, hang, pai, lin):
    """Check the supplied clique, then exclude every one-vertex extension."""
    if len(hang) != jie or len(set(hang)) != jie:
        raise ValueError(f"the witness of order {jie} has the wrong size")

    weizhi = {x: i for i, x in enumerate(pai)}
    dian = set()
    for x in hang:
        if x not in weizhi:
            raise ValueError(f"the witness of order {jie} contains a non-permutation")
        dian.add(weizhi[x])

    for a, b in itertools.combinations(dian, 2):
        if b not in lin[a]:
            raise ValueError(f"the witness of order {jie} is not a clique")

    for x in set(range(120)) - dian:
        keyi = True
        for y in dian:
            if x not in lin[y]:
                keyi = False
                break
        if keyi:
            raise ValueError(f"vertex {x} extends the clique of order {jie}")


def main():
    fen, tuan = dushuju()
    chafen(fen)
    pai, lin, bianshu = jiantu(fen)
    if len(pai) != 120:
        raise ValueError("S_5 must contain 120 permutations")
    for jie in xuyao:
        chatuan(jie, tuan[jie], pai, lin)

    shunxu = ", ".join(map(str, xuyao))
    print(f"completion graph: 120 vertices, {bianshu} edges")
    print(f"verified inclusion-maximal clique orders: {shunxu}")
    print("PASS")


if __name__ == "__main__":
    main()

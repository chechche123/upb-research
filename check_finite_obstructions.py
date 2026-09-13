
from itertools import combinations, permutations, product


M0 = ((0, 1), (2, 3), (4, 5))
HANG = tuple(range(13))
BIAN = tuple((i, j) for i in HANG for j in range(i + 1, 13))
BIANHAO = {e: i for i, e in enumerate(BIAN)}


def yao(tiaojian, hua):
    if not tiaojian:
        raise AssertionError(hua)


def wanquanpeidui(dian):
    if not dian:
        yield ()
        return

    a = dian[0]
    for b in dian[1:]:
        sheng = tuple(x for x in dian if x not in (a, b))
        for m in wanquanpeidui(sheng):
            e = (min(a, b), max(a, b))
            yield tuple(sorted((e,) + m))


def shengchengpeidui(diangeshu):
    """Enumerate three-edge matchings, with consecutive labels for new vertices."""
    for xindian in range(7):
        xin = tuple(range(diangeshu, diangeshu + xindian))
        for jiu in combinations(range(diangeshu), 6 - xindian):
            yield from wanquanpeidui(jiu + xin)


def youhengjie3(m1, m2, m3):
    """Test whether one edge from each matching can be chosen disjointly."""
    for e1 in m1:
        yong = set(e1)
        for e2 in m2:
            if not yong.isdisjoint(e2):
                continue
            yong2 = yong | set(e2)
            if any(yong2.isdisjoint(e3) for e3 in m3):
                return True
    return False


def keyijiaru(zu, m):
    return all(
        not youhengjie3(zu[i], zu[j], m)
        for i, j in combinations(range(len(zu)), 2)
    )


def M0zidong():
    daan = []
    for paixu in permutations(range(3)):
        for fan in product(range(2), repeat=3):
            f = {}
            for i, (a, b) in enumerate(M0):
                c, d = M0[paixu[i]]
                if fan[i]:
                    c, d = d, c
                f[a] = c
                f[b] = d
            daan.append(f)
    return daan


M0ZIDONG = M0zidong()


def guiyixinbian(edges, qidian=6):
    edges = tuple(sorted(tuple(sorted(e)) for e in edges))
    biao = {}
    xiayige = qidian
    daan = []
    for e in edges:
        ne = []
        for x in e:
            if x >= qidian:
                if x not in biao:
                    biao[x] = xiayige
                    xiayige += 1
                x = biao[x]
            ne.append(x)
        daan.append(tuple(sorted(ne)))
    return tuple(sorted(daan))


def dierbiaozhun(m):
    """Normalize the second matching under the 48 automorphisms of M0."""
    houxuan = []
    for f in M0ZIDONG:
        ne = [(f.get(a, a), f.get(b, b)) for a, b in m]
        houxuan.append(guiyixinbian(ne))
    return min(houxuan)


def soupeiduizu():
    diers = sorted({dierbiaozhun(m) for m in shengchengpeidui(6)})
    yao(len(diers) == 27, "the second level must have 27 lists")

    cengshu = [0] * 6
    sizu = []

    def jia(zu):
        n = len(zu)
        cengshu[n] += 1
        if n == 4:
            sizu.append(tuple(zu))
        if n == 5:
            raise AssertionError("a family of five matchings without a transversal exists")

        diangeshu = 1 + max(x for m in zu for e in m for x in e)
        for m in shengchengpeidui(diangeshu):
            if keyijiaru(zu, m):
                jia(zu + [m])

    # Only the second matching is reduced by symmetry. Later levels are exhaustive.
    for m in diers:
        jia([M0, m])

    shiji = tuple(cengshu[2:6])
    yao(shiji == (27, 62, 90, 0), "unexpected matching list counts")
    yao(len(sizu) == 90, "the fourth level must contain 90 lists")
    zuidadian = max(x for zu in sizu for m in zu for e in m for x in e)
    yao(zuidadian < 13, "a matching list uses too many vertices")
    return shiji, tuple(sizu)


def zhuruyouhengjie3(a, b, c):
    for i, j, k in permutations(range(3)):
        wai = (a[i], b[j], c[k])
        if len(set(wai)) == 3:
            return True
    return False


def shengchengzhuru(diangeshu):
    jiudian = range(3, diangeshu)
    for xindian in range(4):
        xin = tuple(range(diangeshu, diangeshu + xindian))
        for jiu in combinations(jiudian, 3 - xindian):
            yield from permutations(jiu + xin)


def zhurukeyijiaru(zu, f):
    return all(
        not zhuruyouhengjie3(zu[i], zu[j], f)
        for i, j in combinations(range(len(zu)), 2)
    )


def shi6huan(zu):
    cishu = {f: zu.count(f) for f in set(zu)}
    if sorted(cishu.values()) != [2, 2]:
        return False
    a, b = cishu
    if set(a) != set(b):
        return False
    if len(set(a)) != 3:
        return False
    return all(a[i] != b[i] for i in range(3))


def souguding3():
    cengshu = [0] * 5

    def jia(zu):
        n = len(zu)
        cengshu[n] += 1
        if n == 4:
            yao(shi6huan(zu), "a fourth-level list lacks the required six-cycle form")
            return

        diangeshu = 1 + max(x for f in zu for x in f)
        for f in shengchengzhuru(diangeshu):
            if zhurukeyijiaru(zu, f):
                jia(zu + [f])

    jia([(3, 4, 5)])
    shiji = tuple(cengshu[1:5])
    yao(shiji == (1, 48, 6, 6), "unexpected fixed-three-set list counts")
    return shiji


def bian(a, b):
    if a < b:
        return a, b
    return b, a


def weizhi(mask):
    while mask:
        zuihou = mask & -mask
        yield zuihou.bit_length() - 1
        mask ^= zuihou


def O3xuanxiang(m):
    """Return the 315 capacity-nine columns for an O_3 matching.

    Each pair contains the covered-edge mask and the degree-two row mask.
    """
    duandian = {x for e in m for x in e}
    yao(len(duandian) == 6, "an O_3 matching does not have six endpoints")
    daan = []

    for ludui in range(3):
        huandui = [i for i in range(3) if i != ludui]
        zuo = m[huandui[0]]
        you = m[huandui[1]]
        huanbian = {bian(a, b) for a in zuo for b in you}
        for zhong in set(HANG) - duandian:
            lubian = {bian(x, zhong) for x in m[ludui]}
            dan = tuple(sorted(set(HANG) - duandian - {zhong}))
            for danm in wanquanpeidui(dan):
                qiangzhi = huanbian | lubian | set(danm)
                jiayi = set(zuo) | set(you) | {zhong}
                yao(len(qiangzhi) == 9, "an O_3 column has the wrong capacity")
                yao(len(jiayi) == 5, "an O_3 option has the wrong bonus set")
                bianmask = sum(1 << BIANHAO[e] for e in qiangzhi)
                hangmask = sum(1 << x for x in jiayi)
                daan.append((bianmask, hangmask))

    yao(len(daan) == 315, "an O_3 matching must have 315 options")
    return tuple(daan)


def bujiaobian(xuan1, xuan2):
    daan = []
    for bian1, _ in xuan1:
        mask = 0
        for i, (bian2, _) in enumerate(xuan2):
            if not bian1 & bian2:
                mask |= 1 << i
        daan.append(mask)
    return tuple(daan)


def sou4O3(zu, xuanku, rongku):
    xuan = [xuanku[m] for m in zu]
    rong = {}
    for i in range(4):
        for j in range(i + 1, 4):
            key = (zu[i], zu[j])
            if key not in rongku:
                rongku[key] = bujiaobian(xuan[i], xuan[j])
            rong[i, j] = rongku[key]

    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0

    for i0, (_, z0) in enumerate(xuan[0]):
        c1 += 1
        h1 = rong[0, 1][i0]
        h20 = rong[0, 2][i0]
        h30 = rong[0, 3][i0]

        for i1 in weizhi(h1):
            c2 += 1
            _, z1 = xuan[1][i1]
            h2 = h20 & rong[1, 2][i1]
            h3di = h30 & rong[1, 3][i1]
            liangci2 = z0 & z1

            for i2 in weizhi(h2):
                c3 += 1
                _, z2 = xuan[2][i2]
                # A row cannot have degree two in three of the chosen columns.
                if liangci2 & z2:
                    continue
                liangci3 = liangci2 | (z0 & z2) | (z1 & z2)
                h3 = h3di & rong[2, 3][i2]
                c4 += h3.bit_count()
                for i3 in weizhi(h3):
                    _, z3 = xuan[3][i3]
                    if not liangci3 & z3:
                        raise AssertionError("four compatible columns of capacity 9 exist")

    return c1, c2, c3, c4


def jiancha4O3(sizu):
    suoyoum = {m for zu in sizu for m in zu}
    xuanku = {m: O3xuanxiang(m) for m in suoyoum}
    rongku = {}
    zong = [0, 0, 0, 0]
    for zu in sizu:
        zhezu = sou4O3(zu, xuanku, rongku)
        zong = [a + b for a, b in zip(zong, zhezu)]

    shiji = tuple(zong)
    yao(shiji == (28_350, 2_400_840, 13_296_960, 0), "unexpected O_3 prefix counts")
    return shiji


def bianbuxiangjiao(e1, e2):
    return e1[0] not in e2 and e1[1] not in e2


def lieyouhengjie3(lie):
    for i, j, k in combinations(range(len(lie)), 3):
        for e1, e2, e3 in product(lie[i], lie[j], lie[k]):
            if (
                bianbuxiangjiao(e1, e2)
                and bianbuxiangjiao(e1, e3)
                and bianbuxiangjiao(e2, e3)
            ):
                return True
    return False


def wanquanbian(zuo, you):
    return frozenset(bian(a, b) for a in zuo for b in you)


def chongfuliang(kuai):
    """Count repeated edges: sum of sizes minus the size of the union."""
    zongshu = sum(len(x) for x in kuai)
    butong = set().union(*kuai) if kuai else set()
    return zongshu - len(butong)


def jiancha2S1():
    x, y, a, b = 0, 1, 2, 3
    wai = tuple(range(4, 13))
    san1 = (x, y, a)
    san2 = (x, y, b)
    duixuan = tuple(combinations(range(3), 2))

    zongshu = 0
    wuhengjie = 0
    O3shaoyi = {0: 10**9, 1: 10**9, 2: 10**9}
    S1shaoyi = 10**9

    for p1 in wai:
        S1dui1 = bian(b, p1)
        for p2 in wai:
            S1dui2 = bian(a, p2)
            S1kuai1 = wanquanbian(san1, S1dui1)
            S1kuai2 = wanquanbian(san2, S1dui2)

            for q1, r1 in permutations(wai, 2):
                O3m1 = (bian(a, b), bian(x, q1), bian(y, r1))
                for q2, r2 in permutations(wai, 2):
                    zongshu += 1
                    O3m2 = (bian(a, b), bian(x, q2), bian(y, r2))
                    lies = ((S1dui1,), (S1dui2,), O3m1, O3m2)
                    if lieyouhengjie3(lies):
                        continue

                    wuhengjie += 1
                    O3ms = (O3m1, O3m2)
                    # A column below capacity contributes no forced K_{2,2}.
                    for shao in range(3):
                        man = 2 - shao
                        for qiyong in combinations(range(2), man):
                            for xuanzes in product(range(3), repeat=man):
                                kuai = [S1kuai1, S1kuai2]
                                for oi, xuan in zip(qiyong, xuanzes):
                                    u, v = duixuan[xuan]
                                    kuai.append(wanquanbian(O3ms[oi][u], O3ms[oi][v]))
                                O3shaoyi[shao] = min(O3shaoyi[shao], chongfuliang(kuai))

                    for xuan1, xuan2 in product(range(3), repeat=2):
                        u1, v1 = duixuan[xuan1]
                        u2, v2 = duixuan[xuan2]
                        kuai = [
                            S1kuai2,
                            wanquanbian(O3m1[u1], O3m1[v1]),
                            wanquanbian(O3m2[u2], O3m2[v2]),
                        ]
                        S1shaoyi = min(S1shaoyi, chongfuliang(kuai))

    yao(zongshu == 419_904, "unexpected labelled support count")
    yao(wuhengjie == 4_104, "unexpected count without a transversal")
    sunshi = (O3shaoyi[0], O3shaoyi[1], O3shaoyi[2], S1shaoyi)
    yao(sunshi == (5, 3, 1, 1), "unexpected two-S_1 losses")
    return zongshu, wuhengjie, sunshi


def jiancha1S1shaoyi():
    m1 = (bian(0, 3), bian(1, 4), bian(2, 5))
    m2 = (bian(0, 4), bian(1, 5), bian(2, 3))
    O3ms = (m1, m1, m2, m2)
    duixuan = tuple(combinations(range(3), 2))

    geshu = []
    zuixiao = []
    for shao in (0, 1):
        man = 4 - shao
        c = 0
        xiao = 10**9
        for qiyong in combinations(range(4), man):
            for xuanzes in product(range(3), repeat=man):
                c += 1
                kuai = []
                for oi, xuan in zip(qiyong, xuanzes):
                    u, v = duixuan[xuan]
                    kuai.append(wanquanbian(O3ms[oi][u], O3ms[oi][v]))
                xiao = min(xiao, chongfuliang(kuai))
        geshu.append(c)
        zuixiao.append(xiao)

    geshu = tuple(geshu)
    zuixiao = tuple(zuixiao)
    yao(geshu == (81, 108), "unexpected one-S_1 counts")
    yao(zuixiao == (4, 2), "unexpected one-S_1 losses")
    return geshu, zuixiao


def main():
    peidui_ceng, sizu = soupeiduizu()
    guding_ceng = souguding3()
    O3_ceng = jiancha4O3(sizu)
    biaoh, wuhengjie, liangS1 = jiancha2S1()
    yiS1geshu, yiS1sunshi = jiancha1S1shaoyi()

    print("matching-family prefixes without a transversal:", *peidui_ceng)
    print("fixed three element set prefixes:", *guding_ceng)
    print("four O_3 column search levels:", *O3_ceng)
    print(f"two S_1 labelled supports: {biaoh} -> {wuhengjie}")
    print("two S_1 minimum repetition losses:", *liangS1)
    print(
        "one S_1 below capacity, patterns and minima:",
        f"{yiS1geshu[0]} -> {yiS1sunshi[0]},",
        f"{yiS1geshu[1]} -> {yiS1sunshi[1]}",
    )
    print("PASS")


if __name__ == "__main__":
    main()

import re
from functools import lru_cache
from pathlib import Path


LUJING = Path(__file__).resolve().parent
SHUJU_WENJIAN = LUJING / "finite_uoms.txt"
XUYAO = (
    (17, 12), (18, 12), (19, 12), (21, 12), (22, 12), (23, 12),
    (18, 14), (19, 14), (21, 14), (22, 14), (23, 14),
    (21, 16), (22, 16), (23, 16), (25, 16),
    (22, 18), (23, 18), (25, 18), (26, 18), (27, 18), (29, 18),
    (30, 18), (31, 18),
)


def banlv(x):
    if x <= 0:
        raise ValueError("matrix symbols must be positive integers")
    if x % 2 == 1:
        return x + 1
    return x - 1


def xingzhuang(juzhen, mingzi):
    if not juzhen or not juzhen[0]:
        raise ValueError(f"{mingzi}: empty matrix")
    lieshu = len(juzhen[0])
    for hang in juzhen:
        if len(hang) != lieshu:
            raise ValueError(f"{mingzi}: nonrectangular matrix")
        if any(type(x) is not int or x <= 0 for x in hang):
            raise ValueError(f"{mingzi}: matrix symbols must be positive integers")
    return len(juzhen), lieshu


def dushuju(wenjian=SHUJU_WENJIAN):
    hangmen = [
        (i, hang.strip())
        for i, hang in enumerate(Path(wenjian).read_text(encoding="utf-8").splitlines(), 1)
        if hang.strip()
    ]
    guize = "MATE RULE: In each column, 1 <-> 2, 3 <-> 4, 5 <-> 6, and so on."
    if not hangmen or hangmen[0][1] != guize:
        raise ValueError("missing or unexpected mate rule")

    quanbu = []
    weizhi = 1
    if len(hangmen) > 1 and hangmen[1][1].startswith("EXPLICIT MATRICES BELOW:"):
        weizhi += 1
    while weizhi < len(hangmen):
        hanghao, neirong = hangmen[weizhi]
        pipei = re.fullmatch(r"SIZE: ([1-9][0-9]*) x ([1-9][0-9]*)", neirong)
        if pipei is None:
            raise ValueError(f"line {hanghao}: expected a SIZE declaration")
        hangshu, lieshu = map(int, pipei.groups())
        mingzi = f"{hangshu} x {lieshu}"
        if weizhi + 2 >= len(hangmen):
            raise ValueError(f"{mingzi}: incomplete matrix record")
        fangfa = hangmen[weizhi + 1][1]
        if not fangfa.startswith("METHOD:") or not fangfa[7:].strip():
            raise ValueError(f"{mingzi}: missing METHOD declaration")
        if hangmen[weizhi + 2][1] != "MATRIX:":
            raise ValueError(f"{mingzi}: missing MATRIX declaration")
        weizhi += 3
        if len(hangmen) - weizhi < hangshu:
            raise ValueError(f"{mingzi}: fewer rows than declared")
        juzhen = []
        for hanghao, neirong in hangmen[weizhi:weizhi + hangshu]:
            if re.fullmatch(r"[1-9][0-9]*(?:\s+[1-9][0-9]*)*", neirong) is None:
                raise ValueError(f"line {hanghao}: expected positive integer entries")
            juzhen.append(tuple(map(int, neirong.split())))
        if xingzhuang(juzhen, mingzi) != (hangshu, lieshu):
            raise ValueError(f"{mingzi}: matrix size differs from its declaration")
        quanbu.append(tuple(juzhen))
        weizhi += hangshu

    chichun = sorted(xingzhuang(x, "matrix") for x in quanbu)
    if chichun != sorted(XUYAO):
        raise ValueError("expected exactly the 23 matrix sizes listed in the paper")
    return quanbu


def jiancha_zhengjiao(juzhen, mingzi):
    hangshu, lieshu = xingzhuang(juzhen, mingzi)
    for i in range(hangshu):
        for j in range(i + 1, hangshu):
            if not any(banlv(juzhen[i][k]) == juzhen[j][k] for k in range(lieshu)):
                raise AssertionError(
                    f"{mingzi}: rows {i + 1} and {j + 1} are not orthogonal"
                )


def dengjia_xianwei(juzhen):
    _, lieshu = xingzhuang(juzhen, "matrix")
    quanbu = []
    for j in range(lieshu):
        zidian = {}
        for i, hang in enumerate(juzhen):
            x = hang[j]
            zidian[x] = zidian.get(x, 0) | (1 << i)
        quanbu.append(zidian)
    return quanbu


def zhao_fugai(juzhen):
    hangshu, lieshu = xingzhuang(juzhen, "matrix")
    xianwei = dengjia_xianwei(juzhen)

    @lru_cache(maxsize=None)
    def sousuo(weifu, liema):
        if weifu == 0:
            return ()
        if liema == 0:
            return None
        kelie = [j for j in range(lieshu) if (liema >> j) & 1]
        weihang = [i for i in range(hangshu) if (weifu >> i) & 1]
        if len(weihang) <= len(kelie):
            return tuple((j, juzhen[i][j]) for i, j in zip(weihang, kelie))
        shangjie = sum(
            max((kuai & weifu).bit_count() for kuai in xianwei[j].values())
            for j in kelie
        )
        if shangjie < len(weihang):
            return None

        hang = min(
            weihang,
            key=lambda i: sum(
                (xianwei[j][juzhen[i][j]] & weifu).bit_count() for j in kelie
            ),
        )
        shunxu = sorted(
            kelie,
            key=lambda j: (xianwei[j][juzhen[hang][j]] & weifu).bit_count(),
            reverse=True,
        )
        for j in shunxu:
            x = juzhen[hang][j]
            kuai = xianwei[j][x]
            jieguo = sousuo(weifu & ~kuai, liema & ~(1 << j))
            if jieguo is not None:
                return ((j, x),) + jieguo
        return None

    jieguo = sousuo((1 << hangshu) - 1, (1 << lieshu) - 1)
    sousuo.cache_clear()
    return jieguo


def jiancha_uom(juzhen, mingzi):
    jiancha_zhengjiao(juzhen, mingzi)
    fugai = zhao_fugai(juzhen)
    if fugai is not None:
        zhengju = ", ".join(f"column {j + 1}: symbol {x}" for j, x in fugai)
        raise AssertionError(f"{mingzi}: extendible; covering fibres: {zhengju}")


def zhu_cheng_xu():
    quanbu = dushuju()
    for i, juzhen in enumerate(quanbu, 1):
        hangshu, lieshu = xingzhuang(juzhen, "matrix")
        mingzi = f"{hangshu} x {lieshu}"
        jiancha_uom(juzhen, mingzi)
        print(f"{i:02d}/23: {mingzi}: orthogonal and unextendible: PASS", flush=True)
    print("verified all 23 matrices in finite_uoms.txt")
    print("PASS")


if __name__ == "__main__":
    zhu_cheng_xu()

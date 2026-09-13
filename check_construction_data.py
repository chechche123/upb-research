
import itertools
import json
from pathlib import Path


LUJING = Path(__file__).resolve().parent
SHUJU_WENJIAN = LUJING / "construction_data.json"
# A cover uses at most one fibre per column. This value marks no cover.
WUQIONG = 10**6

ERKUAI_HANG = {
    4: [6, 7, 8, 9, 10],
    5: [6, 8, 9, 10],
}

SIKUAI_HANG = {
    4: [8, 9, 10, 12],
    5: [8, 9, 10, 11, 12],
    6: [8, 9, 11, 12],
    7: [8, 10, 11, 12, 13],
}

SIKUAI_CANSU = {
    4: [1, 3, 4],
    5: [2, 3, 4],
    6: [2, 3, 5],
    7: [2, 4, 6],
}

LIUKUAI_HANG = [6, 9, 10, 11, 12, 13, 14, 15, 16]
LIUKUAI_CANSU = [1, 2, 3, 4, 5]
LIANGHAO_ZUIDI = [1, 1, 2, 3, 4, 5, 6, 7, 8, None]
LIANGHAO_BUHAN8 = [1, 2, 3, 4, 5, 6, 7, None, None, None]


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
    return len(juzhen), lieshu


def jiancha_zhengjiao(juzhen, mingzi):
    hangshu, lieshu = xingzhuang(juzhen, mingzi)
    for i in range(hangshu):
        for j in range(i + 1, hangshu):
            zhengjiao = False
            for k in range(lieshu):
                if banlv(juzhen[i][k]) == juzhen[j][k]:
                    zhengjiao = True
                    break
            if not zhengjiao:
                raise AssertionError(
                    f"{mingzi}: rows {i + 1} and {j + 1} are not orthogonal"
                )


def dengjia_xianwei(juzhen):
    """Return the equal-symbol row sets in each column, as bit masks."""
    _, lieshu = xingzhuang(juzhen, "matrix")
    quanbu = []
    for j in range(lieshu):
        zidian = {}
        for i, hang in enumerate(juzhen):
            x = hang[j]
            zidian[x] = zidian.get(x, 0) | (1 << i)
        quanbu.append(sorted(set(zidian.values())))
    return quanbu


def zui_xiao_fugai(juzhen, yunxu=None):
    """For each row set, find its minimum cover using the allowed columns.

    A cover chooses at most one equal-symbol fibre from each column.
    Sets with no cover retain the value WUQIONG.
    """
    hangshu, lieshu = xingzhuang(juzhen, "matrix")
    if yunxu is None:
        yunxu = tuple(range(lieshu))
    if len(set(yunxu)) != len(yunxu):
        raise ValueError("an allowed column was repeated")
    for j in yunxu:
        if j < 0 or j >= lieshu:
            raise ValueError("an allowed column lies outside the matrix")

    daijia = [WUQIONG] * (1 << hangshu)
    daijia[0] = 0
    xianwei = dengjia_xianwei(juzhen)
    for j in yunxu:
        # Read the previous table so column j cannot be used twice.
        jiude = daijia[:]
        for mubiao in range(1 << hangshu):
            zuihao = jiude[mubiao]
            for kuai in xianwei[j]:
                zuihao = min(zuihao, jiude[mubiao & ~kuai] + 1)
            daijia[mubiao] = zuihao
    return daijia


def hang_zhuan_ma(hangmen, hangshu, mingzi):
    if not hangmen or len(hangmen) != len(set(hangmen)):
        raise ValueError(f"{mingzi}: a decomposition block is empty or repeats a row")
    for i in hangmen:
        if i < 1 or i > hangshu:
            raise ValueError(f"{mingzi}: a decomposition block contains an invalid row")
    ma = 0
    for i in hangmen:
        ma |= 1 << (i - 1)
    return ma


def jiancha_fenkuai(fenkuai, hangshu, mingzi):
    bing = 0
    for kuai in fenkuai:
        if kuai <= 0:
            raise ValueError(f"{mingzi}: decomposition blocks must be nonempty")
        if bing & kuai:
            raise ValueError(f"{mingzi}: decomposition blocks overlap")
        bing |= kuai
    if bing != (1 << hangshu) - 1:
        raise ValueError(f"{mingzi}: decomposition blocks do not partition the rows")


def fenhua_canshu(daijia, fenkuai):
    """Take the minimum cover number over unions of r blocks, for each r."""
    geshu = len(fenkuai)
    jieguo = []
    for r in range(1, geshu):
        zuihao = WUQIONG
        for zuhe in itertools.combinations(range(geshu), r):
            mubiao = 0
            for i in zuhe:
                mubiao |= fenkuai[i]
            zuihao = min(zuihao, daijia[mubiao])
        jieguo.append(zuihao)
    return jieguo


def neng_fugai(xianwei, mubiao, liema):
    nengdao = {0}
    for j, zhe_lie in enumerate(xianwei):
        if not (liema >> j) & 1:
            continue
        kequ = {kuai & mubiao for kuai in zhe_lie}
        jiude = tuple(nengdao)
        for a in jiude:
            for b in kequ:
                nengdao.add(a | b)
        if mubiao in nengdao:
            return True
    return mubiao in nengdao


def jixiao_zhicheng(juzhen, mubiao):
    """Find the inclusion-minimal column sets that cover the target rows."""
    xianwei = dengjia_xianwei(juzhen)
    jieguo = []
    houxuan = sorted(range(1 << len(xianwei)), key=int.bit_count)
    for liema in houxuan:
        beibaohan = False
        for jiude in jieguo:
            if jiude & liema == jiude:
                beibaohan = True
                break
        if beibaohan:
            continue
        if neng_fugai(xianwei, mubiao, liema):
            jieguo.append(liema)
    return jieguo


def zhicheng_bu(juzhen, fenkuai):
    geshu = len(fenkuai)
    jieguo = set()
    for xuanze in range(1, (1 << geshu) - 1):
        mubiao = 0
        for i, kuai in enumerate(fenkuai):
            if (xuanze >> i) & 1:
                mubiao |= kuai
        for liema in jixiao_zhicheng(juzhen, mubiao):
            jieguo.add((xuanze.bit_count(), liema))
    return jieguo


def chongfu_rongliang(jilu, kuaishu, lieshu):
    """Compute the repetition capacity from disjoint column supports."""
    bu = set()
    for x in jilu:
        bu.update(zhicheng_bu(x["matrix"], x["decomposition_block_masks"]))

    # The table records the largest covered block count for each used support.
    zuiyuan = [-1] * (1 << lieshu)
    zuiyuan[0] = 0
    # Every nonempty support increases the mask size, so this order is enough.
    for yongguo in sorted(range(1 << lieshu), key=int.bit_count):
        if zuiyuan[yongguo] < 0:
            continue
        for gaizhu, liema in bu:
            if yongguo & liema:
                continue
            xinma = yongguo | liema
            zuiyuan[xinma] = max(
                zuiyuan[xinma], min(kuaishu, zuiyuan[yongguo] + gaizhu)
            )
    return max(zuiyuan)


def jiancha_muti(jilu, yingyou_kuaishu):
    mingzi = str(jilu["id"])
    juzhen = jilu["matrix"]
    hangshu, lieshu = xingzhuang(juzhen, mingzi)
    fenkuai = jilu["decomposition_block_masks"]
    if len(fenkuai) != yingyou_kuaishu:
        raise ValueError(f"{mingzi}: wrong number of decomposition blocks")
    jiancha_fenkuai(fenkuai, hangshu, mingzi)
    jiancha_zhengjiao(juzhen, mingzi)
    daijia = zui_xiao_fugai(juzhen)
    if daijia[(1 << hangshu) - 1] <= lieshu:
        raise AssertionError(f"{mingzi}: the input matrix is extendible")
    return daijia, fenhua_canshu(daijia, fenkuai)


def jiancha_erkuai(shuju):
    fenlei = {}
    for x in shuju["two_block_inputs"]:
        mingzi = str(x["id"])
        juzhen = x["matrix"]
        hangshu, lieshu = xingzhuang(juzhen, mingzi)
        fenkuai = []
        for hangmen in x["decomposition_blocks"]:
            fenkuai.append(hang_zhuan_ma(hangmen, hangshu, mingzi))
        y = {
            "id": mingzi,
            "matrix": juzhen,
            "decomposition_block_masks": fenkuai,
        }
        daijia, _ = jiancha_muti(y, 2)
        for kuai in fenkuai:
            if 2 * daijia[kuai] <= lieshu:
                raise AssertionError(f"{mingzi}: the two block condition fails")
        fenlei.setdefault(lieshu, []).append(hangshu)

    if set(fenlei) != {4, 5}:
        raise AssertionError(
            f"two block families: expected columns [4, 5], computed {sorted(fenlei)}"
        )
    for lieshu, yingyou in ERKUAI_HANG.items():
        shiji = sorted(fenlei[lieshu])
        if shiji != yingyou:
            raise AssertionError(
                f"two block N={lieshu}: expected rows {yingyou}, computed {shiji}"
            )
    return fenlei


def jiancha_sikuai(shuju):
    jiating = shuju["four_block_families"]
    lieshu_men = [x["columns"] for x in jiating]
    if len(lieshu_men) != len(set(lieshu_men)) or set(lieshu_men) != {4, 5, 6, 7}:
        raise AssertionError(
            "four block families: columns must occur once each for N=4,5,6,7"
        )

    rongliang = {}
    for jia in jiating:
        lieshu = jia["columns"]
        jilu = jia["inputs"]
        hangmen = []
        for x in jilu:
            hangshu, shiji_lieshu = xingzhuang(x["matrix"], str(x["id"]))
            if shiji_lieshu != lieshu:
                raise ValueError(f"{x['id']}: expected {lieshu} columns")
            jiancha_muti(x, 4)
            hangmen.append(hangshu)

        if sorted(hangmen) != SIKUAI_HANG[lieshu]:
            raise AssertionError(
                f"four block N={lieshu}: the row number list is incomplete"
            )
        zhi = chongfu_rongliang(jilu, 4, lieshu)
        if zhi != 3:
            raise AssertionError(
                f"four block N={lieshu}: repetition capacity is {zhi}, not 3"
            )

        canshu = [WUQIONG, WUQIONG, WUQIONG]
        for x in jilu:
            _, zhege = jiancha_muti(x, 4)
            canshu = [min(a, b) for a, b in zip(canshu, zhege, strict=True)]
            fenkuai = x["decomposition_block_masks"]
            if lieshu == 4:
                for kuai in fenkuai:
                    for liema in jixiao_zhicheng(x["matrix"], kuai):
                        if liema.bit_count() == 1 and not liema & (1 << 3):
                            raise AssertionError(
                                f"{x['id']}: a size 1 support for one block omits column 4"
                            )
            if lieshu == 6:
                for a, b in itertools.combinations(fenkuai, 2):
                    for liema in jixiao_zhicheng(x["matrix"], a | b):
                        if liema.bit_count() == 3 and not liema & (1 << 5):
                            raise AssertionError(
                                f"{x['id']}: a size 3 support for two blocks omits column 6"
                            )
        if canshu != SIKUAI_CANSU[lieshu]:
            raise AssertionError(f"four block N={lieshu}: lower profile is {canshu}")
        rongliang[lieshu] = zhi
    return rongliang


def jiancha_liukuai(shuju):
    jiating = shuju["six_block_family"]
    if jiating["columns"] != 5:
        raise AssertionError("six block family: the column count must be 5")

    canshu_men = []
    hangmen = []
    for x in jiating["inputs"]:
        hangshu, lieshu = xingzhuang(x["matrix"], str(x["id"]))
        if lieshu != 5:
            raise ValueError(f"{x['id']}: wrong column count")
        _, canshu = jiancha_muti(x, 6)
        for a, b in zip(canshu, LIUKUAI_CANSU, strict=True):
            if a < b:
                raise AssertionError(f"{x['id']}: the six block cover condition fails")
        canshu_men.append(canshu)
        hangmen.append(hangshu)

    zuidi = []
    for i in range(5):
        zuidi.append(min(x[i] for x in canshu_men))
    if zuidi != LIUKUAI_CANSU:
        raise AssertionError(f"six block family: lower profile is {zuidi}")
    if sorted(hangmen) != LIUKUAI_HANG:
        raise AssertionError("six block family: the row number list is incomplete")
    return zuidi


def jiancha_banlv_qiquan(juzhen, mingzi):
    for j, zhe_lie in enumerate(zip(*juzhen), start=1):
        fuhao = set(zhe_lie)
        que = sorted(x for x in fuhao if banlv(x) not in fuhao)
        if que:
            raise AssertionError(f"{mingzi}: column {j} omits mates of {que}")


def xianshi_zhi(zhi, lieshu):
    if zhi > lieshu:
        return None
    return zhi


def jiancha_lianghao_lie(shuju):
    jilu = shuju["good_column_input"]
    mingzi = str(jilu["id"])
    juzhen = jilu["matrix"]
    hangshu, lieshu = xingzhuang(juzhen, mingzi)
    if (hangshu, lieshu) != (13, 8):
        raise ValueError(f"{mingzi}: expected a 13 by 8 matrix")
    jiancha_zhengjiao(juzhen, mingzi)
    jiancha_banlv_qiquan(juzhen, mingzi)

    fenkuai = []
    for hangmen in jilu["decomposition_blocks"]:
        fenkuai.append(hang_zhuan_ma(hangmen, hangshu, mingzi))
    jiancha_fenkuai(fenkuai, hangshu, mingzi)
    if len(fenkuai) != 10:
        raise ValueError(f"{mingzi}: expected ten decomposition blocks")
    if jilu["chosen_column"] != 8:
        raise ValueError(f"{mingzi}: the chosen column must be column 8")

    quanlie = tuple(range(lieshu))
    putonglie = tuple(range(7))
    quan_daijia = zui_xiao_fugai(juzhen, quanlie)
    putong_daijia = zui_xiao_fugai(juzhen, putonglie)
    if quan_daijia[(1 << hangshu) - 1] <= lieshu:
        raise AssertionError(f"{mingzi}: the matrix is extendible")

    zuidi_quan = [WUQIONG] * 10
    zuidi_putong = [WUQIONG] * 10
    jishu = 0
    for xuanze in range(1, 1 << 10):
        jishu += 1
        mubiao = 0
        for i, kuai in enumerate(fenkuai):
            if (xuanze >> i) & 1:
                mubiao |= kuai
        k = xuanze.bit_count()
        if quan_daijia[mubiao] < k - 1:
            raise AssertionError(
                f"{mingzi}: the all columns inequality fails at block mask {xuanze}"
            )
        if putong_daijia[mubiao] < k:
            raise AssertionError(
                f"{mingzi}: the inequality without the chosen column fails at block mask {xuanze}"
            )
        zuidi_quan[k - 1] = min(zuidi_quan[k - 1], quan_daijia[mubiao])
        zuidi_putong[k - 1] = min(zuidi_putong[k - 1], putong_daijia[mubiao])

    quan = [xianshi_zhi(x, len(quanlie)) for x in zuidi_quan]
    putong = [xianshi_zhi(x, len(putonglie)) for x in zuidi_putong]
    if quan != LIANGHAO_ZUIDI:
        raise AssertionError(f"{mingzi}: the all columns minima are {quan}")
    if putong != LIANGHAO_BUHAN8:
        raise AssertionError(f"{mingzi}: the minima without the chosen column are {putong}")
    if jishu != 1023:
        raise AssertionError(f"{mingzi}: wrong number of block unions")
    return quan, putong, jishu


def main():
    shuju = json.loads(SHUJU_WENJIAN.read_text(encoding="utf-8"))
    yingyou = {
        "two_block_inputs",
        "four_block_families",
        "six_block_family",
        "good_column_input",
    }
    if set(shuju) != yingyou:
        raise ValueError("construction_data.json has unexpected top-level fields")

    erkuai = jiancha_erkuai(shuju)
    sikuai = jiancha_sikuai(shuju)
    liukuai = jiancha_liukuai(shuju)
    quan, putong, jishu = jiancha_lianghao_lie(shuju)

    erkuai_wenzi = "; ".join(
        f"N={n}: {sorted(hangmen)}" for n, hangmen in sorted(erkuai.items())
    )
    sikuai_wenzi = ", ".join(
        f"N={n}: {zhi}" for n, zhi in sorted(sikuai.items())
    )
    print(f"two block input rows: {erkuai_wenzi}")
    print(f"four block repetition capacities: {sikuai_wenzi}")
    print(f"six block lower profile: {liukuai}")
    print(f"good column X8: {jishu} block unions; minima {quan} / {putong}")
    print("PASS")


if __name__ == "__main__":
    main()

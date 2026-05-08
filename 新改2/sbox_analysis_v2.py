
from itertools import permutations
import time
import random

# ============================================================
# 第一部分：基础工具函数
# ============================================================

def int_to_bits(val, n):
    return [(val >> (n - 1 - i)) & 1 for i in range(n)]

def bits_to_int(bits):
    val = 0
    for b in bits:
        val = (val << 1) | b
    return val

def hamming_weight(val):
    count = 0
    while val:
        count += val & 1
        val >>= 1
    return count

def inner_product(a, b, n):
    return hamming_weight(a & b) % 2

# ============================================================
# 第二部分：S盒拆分与ANF求解
# ============================================================

def sbox_to_component_truth_tables(sbox, n, m):
    size = 1 << n
    components = []
    for bit in range(m):
        tt = [(sbox[x] >> (m - 1 - bit)) & 1 for x in range(size)]
        components.append(tt)
    return components

def mobius_transform(truth_table, n):
    size = 1 << n
    anf = truth_table[:]
    for i in range(n):
        step = 1 << i
        for j in range(size):
            if j & step:
                anf[j] ^= anf[j ^ step]
    return anf

def anf_to_string(anf_coeffs, n, var_names=None):
    if var_names is None:
        var_names = [f"x{i}" for i in range(n)]
    terms = []
    for k in range(1 << n):
        if anf_coeffs[k] == 1:
            if k == 0:
                terms.append("1")
            else:
                terms.append("".join(var_names[i] for i in range(n) if (k >> i) & 1))
    return " ⊕ ".join(terms) if terms else "0"

def compute_algebraic_degree(anf_coeffs, n):
    return max((hamming_weight(k) for k in range(1 << n) if anf_coeffs[k] == 1), default=0)

def count_anf_terms(anf_coeffs):
    return sum(anf_coeffs)

# ============================================================
# 第三部分：Walsh-Hadamard变换与非线性度
# ============================================================

def walsh_hadamard_transform(truth_table, n):
    size = 1 << n
    wht = [1 - 2 * x for x in truth_table]
    for i in range(n):
        step = 1 << i
        for j in range(size):
            if not (j & step):
                u, v = wht[j], wht[j | step]
                wht[j], wht[j | step] = u + v, u - v
    return wht

def nonlinearity_boolean(truth_table, n):
    wht = walsh_hadamard_transform(truth_table, n)
    return (1 << (n - 1)) - max(abs(w) for w in wht) // 2

def nonlinearity_sbox(sbox, n, m):
    size = 1 << n
    min_nl = float('inf')
    for v in range(1, 1 << m):
        tt = [inner_product(v, sbox[x], m) for x in range(size)]
        nl = nonlinearity_boolean(tt, n)
        if nl < min_nl:
            min_nl = nl
    return min_nl

# ============================================================
# 第四部分：差分均匀性
# ============================================================

def difference_distribution_table(sbox, n):
    size = 1 << n
    ddt = [[0] * size for _ in range(size)]
    for dx in range(size):
        for x in range(size):
            ddt[dx][sbox[x] ^ sbox[x ^ dx]] += 1
    return ddt

def differential_uniformity(sbox, n):
    ddt = difference_distribution_table(sbox, n)
    size = 1 << n
    return max(ddt[dx][dy] for dx in range(1, size) for dy in range(size))

# ============================================================
# 第五部分：线性逼近表
# ============================================================

def linear_approximation_table(sbox, n, m):
    size_in, size_out = 1 << n, 1 << m
    lat = [[0] * size_out for _ in range(size_in)]
    for a in range(size_in):
        for b in range(size_out):
            count = sum(1 for x in range(size_in)
                        if inner_product(a, x, n) == inner_product(b, sbox[x], m))
            lat[a][b] = count - (size_in >> 1)
    return lat

def max_linear_approximation(sbox, n, m):
    lat = linear_approximation_table(sbox, n, m)
    return max(abs(lat[a][b]) for a in range(1 << n) for b in range(1, 1 << m))

# ============================================================
# 第六部分：不动点与雪崩特性
# ============================================================

def fixed_points(sbox, n):
    return sum(1 for x in range(1 << n) if sbox[x] == x)

def sac_distance(sbox, n, m):
    size = 1 << n
    components = sbox_to_component_truth_tables(sbox, n, m)
    distances = []
    for tt in components:
        max_dist = max(
            abs(sum(1 for x in range(size) if tt[x] ^ tt[x ^ (1 << (n-1-bit))]) - (size >> 1))
            for bit in range(n)
        )
        distances.append(max_dist)
    return distances

# ============================================================
# 第七部分：隐式代数表示（新增）
# ============================================================

def compute_implicit_representation(sbox, n):
    """
    计算S盒的隐式代数表示。
    对每个分量 yi = fi(x1,...,xn)，移项得到隐式方程 fi(x) ⊕ yi = 0。
    """
    components = sbox_to_component_truth_tables(sbox, n, n)
    x_vars = [f"x{i+1}" for i in range(n)]
    y_vars = [f"y{i+1}" for i in range(n)]
    equations = []

    for i in range(n):
        anf = mobius_transform(components[i], n)
        deg = compute_algebraic_degree(anf, n)

        terms = []
        for k in range(1 << n):
            if anf[k] == 1:
                if k == 0:
                    terms.append("1")
                else:
                    terms.append("".join(x_vars[b] for b in range(n) if (k >> b) & 1))
        terms.append(y_vars[i])  # 移项

        eq_str = " ⊕ ".join(terms) + " = 0"
        is_linear = (deg <= 1)

        equations.append({
            'index': i,
            'anf': anf,
            'degree': deg,
            'is_linear': is_linear,
            'equation_str': eq_str,
        })

    num_linear = sum(1 for eq in equations if eq['is_linear'])
    cls = 'A' if num_linear == 0 else ('C' if num_linear == n else 'B')

    return {
        'equations': equations,
        'num_linear_equations': num_linear,
        'implicit_class': cls,
    }

def print_implicit_representation(sbox, n, name="S-box", file=None):
    result = compute_implicit_representation(sbox, n)
    lines = [f"\n{'='*65}", f"隐式代数表示：{name}", f"{'='*65}"]
    for eq in result['equations']:
        tag = "（线性）" if eq['is_linear'] else f"（次数{eq['degree']}）"
        lines.append(f"  f{eq['index']+1}: {eq['equation_str']}  {tag}")
    lines.append(f"\n  含线性方程数 k = {result['num_linear_equations']}")
    if n == 3:
        lines.append(f"  隐式分类：{result['implicit_class']}类")
    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n")
    return result

# ============================================================
# 第八部分：完整S盒分析
# ============================================================

def analyze_sbox(sbox, n, m, name="S-box"):
    components = sbox_to_component_truth_tables(sbox, n, m)
    var_names = [f"x{i}" for i in range(n)]
    anf_list, anf_strings, degrees, term_counts = [], [], [], []
    for comp in components:
        anf = mobius_transform(comp, n)
        anf_list.append(anf)
        anf_strings.append(anf_to_string(anf, n, var_names))
        degrees.append(compute_algebraic_degree(anf, n))
        term_counts.append(count_anf_terms(anf))
    return {
        'name': name, 'n': n, 'm': m, 'truth_table': sbox,
        'anf_strings': anf_strings, 'degrees': degrees,
        'sbox_degree': max(degrees), 'term_counts': term_counts,
        'nonlinearity': nonlinearity_sbox(sbox, n, m),
        'diff_uniformity': differential_uniformity(sbox, n),
        'max_linear_approx': max_linear_approximation(sbox, n, m),
        'fixed_points': fixed_points(sbox, n),
        'sac_distances': sac_distance(sbox, n, m),
    }

def print_results(results, file=None):
    lines = ["=" * 70, f"S盒名称：{results['name']}",
             f"规模：{results['n']}×{results['m']}",
             f"真值表：{[hex(x) for x in results['truth_table']]}", "-" * 70,
             "【代数正规型（ANF）】"]
    for i in range(results['m']):
        lines.append(f"  y{i} = {results['anf_strings'][i]}")
    lines.append("\n【代数次数】")
    for i in range(results['m']):
        lines.append(f"  deg(y{i}) = {results['degrees'][i]}，项数 = {results['term_counts'][i]}")
    lines += [f"  S盒代数次数 = {results['sbox_degree']}",
              f"\n【非线性度】 NL = {results['nonlinearity']}",
              f"【差分均匀性】 δ = {results['diff_uniformity']}",
              f"【最大线性逼近值】 |LAT_max| = {results['max_linear_approx']}",
              f"【最大线性逼近优势】 A_max = {results['max_linear_approx']}/{1 << results['n']}",
              f"【不动点个数】 {results['fixed_points']}",
              f"【严格雪崩距离】 {results['sac_distances']}", "=" * 70]
    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n\n")

# ============================================================
# 第九部分：3×3双射S盒批量枚举与分类
# ============================================================

def enumerate_3x3_bijective_sboxes():
    print("\n" + "=" * 70)
    print("正在枚举所有3×3双射S盒（共40320个）...")
    print("=" * 70)
    start_time = time.time()
    categories = {}
    total = 0
    for perm in permutations(range(8)):
        sbox = list(perm)
        total += 1
        nl = nonlinearity_sbox(sbox, 3, 3)
        du = differential_uniformity(sbox, 3)
        comps = sbox_to_component_truth_tables(sbox, 3, 3)
        deg = max(compute_algebraic_degree(mobius_transform(c, 3), 3) for c in comps)
        key = (nl, du, deg)
        if key not in categories:
            categories[key] = {'count': 0, 'example': sbox}
        categories[key]['count'] += 1
        if total % 5000 == 0:
            print(f"  已处理 {total}/40320 ...")
    print(f"\n枚举完成，耗时 {time.time()-start_time:.1f} 秒")
    return categories, total

def print_3x3_classification(categories, total, file=None):
    lines = ["\n" + "=" * 70,
             "3×3双射S盒分类结果（按非线性度、差分均匀性、代数次数分类）",
             f"总数：{total}", "=" * 70,
             f"{'非线性度':>8} {'差分均匀性':>10} {'代数次数':>8} {'数量':>8} {'占比':>8} {'示例S盒'}",
             "-" * 70]
    for key in sorted(categories.keys(), reverse=True):
        nl, du, deg = key
        count = categories[key]['count']
        lines.append(f"{nl:>8} {du:>10} {deg:>8} {count:>8} {count/total*100:>7.2f}% {categories[key]['example']}")
    opt = sum(categories[k]['count'] for k in categories if k == (2, 2, 2))
    lines += ["-" * 70, f"最优S盒（NL=2, δ=2, deg=2）数量：{opt}，占比：{opt/total*100:.2f}%", "=" * 70]
    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n\n")

# ============================================================
# 第十部分：任务1 —— 4×4 S盒16个最优等价类隐式表示
# ============================================================

def find_16_optimal_classes(seed=42):
    """
    随机枚举搜索16个仿射等价意义下不同的最优4×4 S盒代表元。
    仿射等价不变量：DDT值分布 + LAT绝对值分布 + 不动点数 + 对合性。
    平均只需约1万次尝试即可找齐16类。
    """
    def invariant_fp(sbox):
        ddt = difference_distribution_table(sbox, 4)
        lat = linear_approximation_table(sbox, 4, 4)
        ddt_v = tuple(sorted(ddt[dx][dy] for dx in range(1,16) for dy in range(16)))
        lat_v  = tuple(sorted(abs(lat[a][b]) for a in range(16) for b in range(1,16)))
        fp  = fixed_points(sbox, 4)
        inv = int(all(sbox[sbox[x]] == x for x in range(16)))
        return (ddt_v, lat_v, fp, inv)

    random.seed(seed)
    perm = list(range(16))
    found = {}
    attempts = 0
    print("正在搜索16个最优等价类代表元（NL=4, δ=4, deg=3）...")
    while len(found) < 16:
        random.shuffle(perm)
        sbox = perm[:]
        attempts += 1
        if nonlinearity_sbox(sbox, 4, 4) != 4:
            continue
        if differential_uniformity(sbox, 4) != 4:
            continue
        comps = sbox_to_component_truth_tables(sbox, 4, 4)
        if max(compute_algebraic_degree(mobius_transform(c,4),4) for c in comps) != 3:
            continue
        key = invariant_fp(sbox)
        if key not in found:
            found[key] = sbox[:]
            print(f"  第{len(found):>2}类（第{attempts}次）: {[hex(x) for x in sbox]}")
    print(f"搜索完成，共尝试 {attempts} 次。\n")
    return list(found.values())


def task1_16class_implicit(file=None):
    """任务1：计算4×4 S盒16个最优等价类代表元的隐式代数表示"""
    lines = ["\n" + "=" * 70,
             "【任务1】4×4 S盒16个最优等价类代表元的隐式代数表示",
             "（NL=4, δ=4, deg=3；代表元由随机枚举搜索得到）",
             "=" * 70]

    reps = find_16_optimal_classes(seed=42)
    summary = []

    for idx, sbox in enumerate(reps, 1):
        nl  = nonlinearity_sbox(sbox, 4, 4)
        du  = differential_uniformity(sbox, 4)
        comps = sbox_to_component_truth_tables(sbox, 4, 4)
        deg = max(compute_algebraic_degree(mobius_transform(c,4),4) for c in comps)
        fp  = fixed_points(sbox, 4)
        inv = all(sbox[sbox[x]] == x for x in range(16))
        impl = compute_implicit_representation(sbox, 4)

        lines.append(f"\n--- 等价类 C{idx} ---")
        lines.append(f"真值表: {[hex(x) for x in sbox]}")
        lines.append(f"指标: NL={nl}, δ={du}, deg={deg}, 不动点={fp}, 对合={'是' if inv else '否'}")
        lines.append("隐式方程组：")
        for eq in impl['equations']:
            tag = "（线性）" if eq['is_linear'] else f"（次数{eq['degree']}）"
            lines.append(f"  f{eq['index']+1}: {eq['equation_str']}  {tag}")
        lines.append(f"含线性方程数 k = {impl['num_linear_equations']}")

        summary.append({'idx': idx, 'nl': nl, 'du': du, 'deg': deg,
                        'fp': fp, 'inv': inv, 'k': impl['num_linear_equations']})

    lines += ["\n\n" + "=" * 70, "16个最优等价类隐式表示汇总表", "=" * 70,
              f"{'类别':>5}  {'NL':>4}  {'δ':>4}  {'deg':>4}  {'不动点':>6}  {'对合':>5}  {'线性方程数k':>10}",
              "-" * 55]
    for r in summary:
        lines.append(f"  C{r['idx']:<3}  {r['nl']:>4}  {r['du']:>4}  {r['deg']:>4}  "
                     f"{r['fp']:>6}  {'是' if r['inv'] else '否':>5}  {r['k']:>10}")

    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n\n")

# ============================================================
# 第十一部分：任务2 —— 3×3双射S盒全部隐式表示枚举分类
# ============================================================

def task2_3x3_implicit_enumeration(file=None):
    """
    枚举全部40320个3×3双射S盒的隐式表示，
    按（隐式类别, NL, δ）统计分布，并给出各组代表例的完整方程。
    """
    print("\n" + "=" * 70)
    print("【任务2】3×3双射S盒全部隐式表示枚举与分类（共40320个）")
    print("=" * 70)
    start = time.time()

    struct_stats = {}
    class_examples = {}
    total = 0

    for perm in permutations(range(8)):
        sbox = list(perm)
        total += 1
        nl = nonlinearity_sbox(sbox, 3, 3)
        du = differential_uniformity(sbox, 3)
        impl = compute_implicit_representation(sbox, 3)
        k = impl['num_linear_equations']
        key = (k, nl, du)
        if key not in struct_stats:
            struct_stats[key] = {'count': 0, 'k': k, 'nl': nl, 'du': du,
                                  'class': impl['implicit_class']}
        struct_stats[key]['count'] += 1
        if key not in class_examples:
            class_examples[key] = (sbox, impl)
        if total % 8000 == 0:
            print(f"  已处理 {total}/40320 ... ({time.time()-start:.1f}s)")

    print(f"\n枚举完成，耗时 {time.time()-start:.1f} 秒")

    lines = ["\n3×3双射S盒隐式表示结构分布表（共40320个）", "=" * 65,
             f"{'隐式类':>6}  {'k(线性方程数)':>13}  {'NL':>4}  {'δ':>4}  {'数量':>8}  {'占比':>8}",
             "-" * 60]
    for key in sorted(struct_stats.keys()):
        k, nl, du = key
        info = struct_stats[key]
        lines.append(f"  {info['class']}类  {k:>13}  {nl:>4}  {du:>4}  "
                     f"{info['count']:>8}  {info['count']/total*100:>7.2f}%")
    lines += ["-" * 60, f"{'合计':>35}  {total:>8}  {'100.00%':>8}",
              "\n\n各组代表S盒的完整隐式方程展示：", "=" * 65]

    for key in sorted(struct_stats.keys()):
        k, nl, du = key
        sbox_ex, impl_ex = class_examples[key]
        info = struct_stats[key]
        lines.append(f"\n--- {info['class']}类（k={k}, NL={nl}, δ={du}，"
                     f"共{info['count']}个，占{info['count']/total*100:.2f}%）---")
        lines.append(f"代表真值表: {sbox_ex}")
        lines.append("隐式方程组：")
        for eq in impl_ex['equations']:
            tag = "（线性）" if eq['is_linear'] else f"（次数{eq['degree']}）"
            lines.append(f"  f{eq['index']+1}: {eq['equation_str']}  {tag}")

    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n\n")
    return struct_stats

# ============================================================
# 第十二部分：主程序
# ============================================================

def main():
    output_file = open("results_v2.txt", "w", encoding="utf-8")
    print("=" * 70)
    print("小规模S盒的代数表示及密码学性质计算（含隐式表示扩展）")
    print("=" * 70)
    output_file.write("小规模S盒的代数表示及密码学性质计算（含隐式表示扩展）\n\n")

    # 一、3×3 S盒分析
    print("\n【一、3×3 S盒分析】\n")
    output_file.write("【一、3×3 S盒分析】\n\n")
    for sbox, name in [([5,3,6,1,0,7,2,4], "S1 (论文实例1)"),
                       ([0,1,3,6,7,4,5,2], "S2 (论文实例2)")]:
        r = analyze_sbox(sbox, 3, 3, name)
        print_results(r, output_file)
        print_implicit_representation(sbox, 3, name, output_file)

    # 二、4×4 S盒分析
    print("\n【二、4×4 S盒分析】\n")
    output_file.write("【二、4×4 S盒分析】\n\n")
    present = [0xC,0x5,0x6,0xB,0x9,0x0,0xA,0xD,0x3,0xE,0xF,0x8,0x4,0x7,0x1,0x2]
    gift    = [0x1,0xA,0x4,0xC,0x6,0xF,0x3,0x9,0x2,0xD,0xB,0x7,0x5,0x0,0x8,0xE]
    midori  = [0xC,0xA,0xD,0x3,0xE,0xB,0xF,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6]
    rp = analyze_sbox(present, 4, 4, "PRESENT S盒")
    rg = analyze_sbox(gift,    4, 4, "GIFT S盒")
    rm = analyze_sbox(midori,  4, 4, "Midori Sb0")
    for r, s in [(rp, present), (rg, gift), (rm, midori)]:
        print_results(r, output_file)
        print_implicit_representation(s, 4, r['name'], output_file)

    # 三、PRESENT差分分布表
    print("\n【三、PRESENT S盒差分分布表】\n")
    output_file.write("【三、PRESENT S盒差分分布表】\n\n")
    ddt = difference_distribution_table(present, 4)
    header = "Δx\\Δy " + " ".join(f"{i:>3X}" for i in range(16))
    print(header); output_file.write(header + "\n")
    for dx in range(16):
        row = f"  {dx:>3X}  " + " ".join(f"{ddt[dx][dy]:>3}" for dy in range(16))
        print(row); output_file.write(row + "\n")
    output_file.write("\n")

    # 四、PRESENT线性逼近表
    print("\n【四、PRESENT S盒线性逼近表】\n")
    output_file.write("【四、PRESENT S盒线性逼近表】\n\n")
    lat = linear_approximation_table(present, 4, 4)
    header = " a\\b  " + " ".join(f"{i:>3X}" for i in range(16))
    print(header); output_file.write(header + "\n")
    for a in range(16):
        row = f"  {a:>3X}  " + " ".join(f"{lat[a][b]:>3}" for b in range(16))
        print(row); output_file.write(row + "\n")
    output_file.write("\n")

    # 五、4×4 S盒综合对比
    print("\n【五、4×4 S盒综合指标对比】\n")
    output_file.write("【五、4×4 S盒综合指标对比】\n\n")
    header = f"{'S盒':>12} {'非线性度':>8} {'差分均匀性':>10} {'代数次数':>8} {'|LAT_max|':>10} {'不动点':>6}"
    print(header); output_file.write(header + "\n")
    print("-"*60); output_file.write("-"*60+"\n")
    for r in [rp, rg, rm]:
        row = (f"{r['name']:>12} {r['nonlinearity']:>8} {r['diff_uniformity']:>10} "
               f"{r['sbox_degree']:>8} {r['max_linear_approx']:>10} {r['fixed_points']:>6}")
        print(row); output_file.write(row + "\n")
    output_file.write("\n")

    # 六、3×3 批量枚举（指标分类）
    print("\n【六、3×3双射S盒批量枚举与分类（指标分类）】")
    output_file.write("【六、3×3双射S盒批量枚举与分类（指标分类）】\n")
    cats, total3 = enumerate_3x3_bijective_sboxes()
    print_3x3_classification(cats, total3, output_file)

    # 七、任务1：16个最优等价类隐式表示
    print("\n【七、4×4 S盒16个最优等价类隐式代数表示】")
    output_file.write("【七、4×4 S盒16个最优等价类隐式代数表示】\n")
    task1_16class_implicit(output_file)

    # 八、任务2：3×3全枚举隐式分类
    print("\n【八、3×3双射S盒全部隐式表示枚举分类】")
    output_file.write("【八、3×3双射S盒全部隐式表示枚举分类】\n")
    task2_3x3_implicit_enumeration(output_file)

    output_file.close()
    print(f"\n✅ 所有结果已保存到 results_v2.txt")

if __name__ == "__main__":
    main()

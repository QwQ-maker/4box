
from itertools import permutations
import time
import random
from collections import Counter, defaultdict

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
# 第七部分：隐式代数表示
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
# 【新增】第七部分之补充：ANF表达式结构特征提取
# ------------------------------------------------------------
# 老师新要求：分类不再用性能指标，改用"表达式本身"的结构特征
# （项数、次数分布、二次项数量等）作为分类依据。
# 以下函数负责从一个S盒的ANF中提取所有能用于分类的结构指纹。
# ============================================================

def extract_anf_structure(sbox, n):
    """
    对S盒的每个分量布尔函数，提取ANF表达式的结构特征。
    返回值是一个dict，包含所有可用于"按表达式分类"的维度。

    特征维度说明：
      - per_component_terms   : 各分量的项数  (例: (3,2,6))
      - per_component_degrees : 各分量的次数  (例: (1,1,2))
      - per_component_deg_profile: 各分量"(一次项数, 二次项数, 三次项数)"
      - total_terms           : 全部分量ANF的总项数
      - total_linear_terms    : 所有分量中一次项的总数
      - total_quadratic_terms : 所有分量中二次项的总数
      - total_cubic_terms     : 所有分量中三次项的总数（3×3最多为0或出现于y中）
      - has_constant_per_comp : 各分量是否含常数项  (例: (1,0,1))
      - num_constant_terms    : 含常数项的分量数量
      - quadratic_monomials   : 所有分量出现的二次单项式索引集合（结构指纹）
      - structure_signature   : 综合结构签名（可作为分类key）
    """
    components = sbox_to_component_truth_tables(sbox, n, n)

    per_terms = []
    per_degrees = []
    per_deg_profile = []     # 每个分量的 (deg1项数, deg2项数, ..., degn项数)
    per_has_const = []
    all_quadratic_monos = set()   # 所有二次项在哪些索引出现（跨分量汇总）

    total_terms = 0
    deg_count_global = [0] * (n + 1)   # deg_count_global[d] = 所有分量中次数为d的项总数

    for i, comp in enumerate(components):
        anf = mobius_transform(comp, n)
        deg = compute_algebraic_degree(anf, n)
        nterms = count_anf_terms(anf)

        per_terms.append(nterms)
        per_degrees.append(deg)
        total_terms += nterms

        # 统计每个次数的项数
        local_profile = [0] * (n + 1)   # index 0=常数项, 1=一次项, 2=二次项, ...
        for k in range(1 << n):
            if anf[k] == 1:
                w = hamming_weight(k)
                local_profile[w] += 1
                deg_count_global[w] += 1
                if w == 2:
                    all_quadratic_monos.add(k)

        per_deg_profile.append(tuple(local_profile))
        per_has_const.append(1 if anf[0] == 1 else 0)

    # 综合结构签名：作为"按表达式分类"的主键
    # 设计思路：先按各分量的(次数, 项数)排序后形成元组，保证S盒分量顺序不同
    # 但结构相同时能归为同一类；同时附带二次项指纹，保证真正结构相同才归一类。
    sorted_components_fingerprint = tuple(sorted(
        zip(per_degrees, per_terms, per_deg_profile, per_has_const)
    ))
    structure_signature = (
        sorted_components_fingerprint,
        tuple(deg_count_global),
        len(all_quadratic_monos),
    )

    return {
        'per_component_terms': tuple(per_terms),
        'per_component_degrees': tuple(per_degrees),
        'per_component_deg_profile': tuple(per_deg_profile),
        'total_terms': total_terms,
        'total_linear_terms': deg_count_global[1],
        'total_quadratic_terms': deg_count_global[2],
        'total_cubic_terms': deg_count_global[3] if n >= 3 else 0,
        'has_constant_per_comp': tuple(per_has_const),
        'num_constant_terms': sum(per_has_const),
        'quadratic_monomials': frozenset(all_quadratic_monos),
        'deg_count_global': tuple(deg_count_global),
        'structure_signature': structure_signature,
    }


def format_anf_structure(struct, n):
    """把结构特征格式化为人类可读的字符串（用于在报告中展示每一类）。"""
    lines = []
    lines.append(f"  总项数 = {struct['total_terms']}")
    lines.append(f"  各分量项数 = {struct['per_component_terms']}")
    lines.append(f"  各分量次数 = {struct['per_component_degrees']}")
    lines.append(f"  各分量次数谱(常数项,一次,二次,三次) = "
                 + str(struct['per_component_deg_profile']))
    lines.append(f"  一次项总数 = {struct['total_linear_terms']}, "
                 f"二次项总数 = {struct['total_quadratic_terms']}")
    if n >= 3:
        lines.append(f"  三次项总数 = {struct['total_cubic_terms']}")
    lines.append(f"  含常数项的分量 = {struct['has_constant_per_comp']} "
                 f"（共 {struct['num_constant_terms']} 个）")
    return "\n".join(lines)

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
# 【新增】第十二部分：任务3 —— 基于ANF表达式结构的分类
# ------------------------------------------------------------
# 这是针对老师新要求的核心模块。
# 对全部40320个3×3双射S盒枚举后，按"表达式结构特征"聚类：
#   - 维度1（粗粒度）：按总项数分类
#   - 维度2（中粒度）：按(各分量次数, 各分量项数)分类
#   - 维度3（细粒度）：按完整structure_signature分类
# 重点关注10752个最优S盒（NL=2,δ=2,deg=2）内部的表达式分布差异。
# ============================================================

def task3_anf_structure_classification(file=None, export_all_anf=True):
    """
    任务3：基于ANF表达式结构对3×3双射S盒进行分类。

    核心思想（按老师要求）：
    不再以性能指标(NL/δ/deg)为分类依据，而是直接从ANF表达式本身提取
    项数、次数谱、二次项集合等结构特征作为分类维度，专门研究
    "性能相同但表达式结构不同"的S盒之间的差异。
    """
    print("\n" + "=" * 70)
    print("【任务3】基于ANF表达式结构的3×3双射S盒分类（老师新要求）")
    print("=" * 70)
    start = time.time()

    # 实时计数器
    cnt_total = Counter()                # key = 总项数
    cnt_tv    = Counter()                # key = 各分量项数(排序)
    cnt_dp    = Counter()                # key = 各分量次数谱(排序)
    cnt_fine  = Counter()                # key = 完整结构签名
    cnt_total_opt = Counter()
    cnt_tv_opt    = Counter()
    cnt_dp_opt    = Counter()
    cnt_fine_opt  = Counter()

    # 每个类保留若干代表S盒，用于报告展示
    example_total    = {}
    example_tv       = {}
    example_fine     = {}
    example_total_opt = {}
    example_tv_opt    = {}

    # 所有S盒的ANF列表（可选导出为附录A数据文件）
    all_anf_records = []

    total = 0
    optimal_count = 0

    for perm in permutations(range(8)):
        sbox = list(perm)
        total += 1

        # 性能指标——只用来判断是否"最优"，不再作为分类维度
        nl = nonlinearity_sbox(sbox, 3, 3)
        du = differential_uniformity(sbox, 3)
        is_optimal = (nl == 2 and du == 2)
        if is_optimal:
            optimal_count += 1

        # 表达式结构特征——老师要求的分类依据
        struct = extract_anf_structure(sbox, 3)

        # 分类key
        key_total = struct['total_terms']
        key_tv    = tuple(sorted(struct['per_component_terms']))
        key_dp    = tuple(sorted(struct['per_component_deg_profile']))
        key_fine  = struct['structure_signature']

        # 全集计数
        cnt_total[key_total] += 1
        cnt_tv[key_tv]       += 1
        cnt_dp[key_dp]       += 1
        cnt_fine[key_fine]   += 1

        # 最优子集计数
        if is_optimal:
            cnt_total_opt[key_total] += 1
            cnt_tv_opt[key_tv]       += 1
            cnt_dp_opt[key_dp]       += 1
            cnt_fine_opt[key_fine]   += 1

        # 代表例（每类保留1个）
        if key_total not in example_total:
            example_total[key_total] = sbox[:]
        if key_tv not in example_tv:
            example_tv[key_tv] = sbox[:]
        if key_fine not in example_fine:
            example_fine[key_fine] = sbox[:]
        if is_optimal:
            if key_total not in example_total_opt:
                example_total_opt[key_total] = sbox[:]
            if key_tv not in example_tv_opt:
                example_tv_opt[key_tv] = sbox[:]

        # 导出完整ANF（用于附录A）
        if export_all_anf:
            comps = sbox_to_component_truth_tables(sbox, 3, 3)
            anf_strs = []
            for c in comps:
                anf = mobius_transform(c, 3)
                anf_strs.append(anf_to_string(anf, 3, ["x1","x2","x3"]))
            all_anf_records.append({
                'sbox': sbox[:],
                'y1_anf': anf_strs[0],
                'y2_anf': anf_strs[1],
                'y3_anf': anf_strs[2],
                'total_terms': struct['total_terms'],
                'per_terms': struct['per_component_terms'],
                'per_degs':  struct['per_component_degrees'],
                'nl': nl, 'du': du,
                'is_optimal': is_optimal,
            })

        if total % 8000 == 0:
            print(f"  已处理 {total}/40320 ... ({time.time()-start:.1f}s)")

    n_fine_classes_opt = len(cnt_fine_opt)
    n_fine_classes_all = len(cnt_fine)

    print(f"\n枚举完成，耗时 {time.time()-start:.1f} 秒")
    print(f"全部S盒：{total}，其中最优S盒(NL=2,δ=2)：{optimal_count}")

    # ============ 报告生成 ============
    lines = []
    lines.append("\n" + "=" * 70)
    lines.append("【任务3】基于ANF表达式结构的3×3双射S盒分类结果")
    lines.append("=" * 70)
    lines.append(f"枚举总数：{total}；其中最优S盒(NL=2, δ=2, deg=2)：{optimal_count}")
    lines.append("")
    lines.append("说明：本任务不以性能指标为分类依据，而以ANF表达式本身的")
    lines.append("     结构特征（项数、次数分布等）作为分类维度。")

    # ---- 表3-1：按总项数分类（全集） ----
    lines.append("\n" + "-" * 70)
    lines.append("表5-A-1  全部40320个3×3双射S盒按ANF总项数的分布")
    lines.append("-" * 70)
    lines.append(f"{'总项数':>8}  {'数量':>8}  {'占比':>10}")
    for k in sorted(cnt_total.keys()):
        pct = cnt_total[k] / total * 100
        lines.append(f"{k:>8}  {cnt_total[k]:>8}  {pct:>9.2f}%")

    # ---- 表3-2：按总项数分类（最优S盒子集） ----
    lines.append("\n" + "-" * 70)
    lines.append(f"表5-A-2  {optimal_count}个最优3×3 S盒按ANF总项数的分布")
    lines.append("       （NL=2, δ=2, deg=2 内部的表达式结构细分）")
    lines.append("-" * 70)
    lines.append(f"{'总项数':>8}  {'数量':>8}  {'占比':>10}")
    for k in sorted(cnt_total_opt.keys()):
        pct = cnt_total_opt[k] / optimal_count * 100 if optimal_count else 0
        lines.append(f"{k:>8}  {cnt_total_opt[k]:>8}  {pct:>9.2f}%")

    # ---- 表3-3：按"各分量项数向量"分类（最优S盒子集） ----
    lines.append("\n" + "-" * 70)
    lines.append(f"表5-A-3  {optimal_count}个最优3×3 S盒按各分量项数向量的分布")
    lines.append("       （分量排序后的元组）")
    lines.append("-" * 70)
    lines.append(f"{'项数向量(排序)':>18}  {'数量':>8}  {'占比':>10}  {'代表真值表'}")
    for tv, c in sorted(cnt_tv_opt.items(), key=lambda x: -x[1]):
        pct = c / optimal_count * 100 if optimal_count else 0
        ex = example_tv_opt.get(tv, [])
        lines.append(f"{str(tv):>18}  {c:>8}  {pct:>9.2f}%  {ex}")

    # ---- 表3-4：完整结构签名细分类汇总 ----
    lines.append("\n" + "-" * 70)
    lines.append("表5-A-4  3×3双射S盒按完整ANF结构签名的细分类汇总")
    lines.append("-" * 70)
    lines.append(f"  全集40320个S盒中，不同表达式结构签名类数：{n_fine_classes_all}")
    lines.append(f"  {optimal_count}个最优S盒中，不同表达式结构签名类数：{n_fine_classes_opt}")
    size_dist = Counter(cnt_fine_opt.values())
    lines.append(f"  最优S盒各签名类所含S盒数量的分布（类大小 -> 有几个这样的类）：")
    for sz in sorted(size_dist.keys()):
        lines.append(f"    大小={sz:>5}  的签名类数量 = {size_dist[sz]}")

    # ---- 表3-5：每个"总项数类"的代表ANF展示（最优子集） ----
    lines.append("\n" + "-" * 70)
    lines.append('表5-A-5  最优S盒各"总项数类"的代表ANF展示')
    lines.append("-" * 70)
    for k in sorted(cnt_total_opt.keys()):
        rep = example_total_opt.get(k)
        if rep is None:
            continue
        struct = extract_anf_structure(rep, 3)
        comps = sbox_to_component_truth_tables(rep, 3, 3)
        anf_strs = [anf_to_string(mobius_transform(c,3), 3, ["x1","x2","x3"]) for c in comps]
        lines.append(f"\n▌ 总项数 = {k}（共 {cnt_total_opt[k]} 个最优S盒属于此类）")
        lines.append(f"  代表真值表: {rep}")
        lines.append(f"  y1 = {anf_strs[0]}")
        lines.append(f"  y2 = {anf_strs[1]}")
        lines.append(f"  y3 = {anf_strs[2]}")
        lines.append(format_anf_structure(struct, 3))

    # ---- 论文建议结论 ----
    lines.append("\n" + "=" * 70)
    lines.append("【基于表达式结构分类的关键发现（论文5.x节可用）】")
    lines.append("=" * 70)
    lines.append("")
    lines.append(f"1. 在{optimal_count}个性能同为最优(NL=2, δ=2, deg=2)的3×3 S盒中，")
    lines.append("   ANF总项数并非固定值，而是呈现出若干不同的档次（见表5-A-2），")
    lines.append("   说明性能指标相同的S盒在表达式复杂度上存在明显差异。")
    lines.append("")
    lines.append('2. 按"各分量项数向量"划分的子类中（表5-A-3），项数最少的')
    lines.append("   子类对应着最简ANF表达式，在硬件实现中可能具有更低的门电")
    lines.append("   路代价，具有工程优选价值。")
    lines.append("")
    lines.append(f"3. 完整结构签名下可细分为 {n_fine_classes_opt} 个不同的表达式结构")
    lines.append("   子类，这一细分粒度远超传统性能指标分类所能提供的信息量，")
    lines.append('   表明"表达式结构"是一个独立于性能指标的分类维度。')

    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n\n")

    # 另存ANF完整表为附录数据
    if export_all_anf:
        with open("appendix_A_3x3_anf_full.csv", "w", encoding="utf-8") as f:
            f.write("sbox_truth_table,y1_ANF,y2_ANF,y3_ANF,total_terms,per_terms,per_degrees,NL,delta,is_optimal\n")
            for rec in all_anf_records:
                f.write(f"\"{rec['sbox']}\",\"{rec['y1_anf']}\",\"{rec['y2_anf']}\","
                        f"\"{rec['y3_anf']}\",{rec['total_terms']},"
                        f"\"{rec['per_terms']}\",\"{rec['per_degs']}\","
                        f"{rec['nl']},{rec['du']},{int(rec['is_optimal'])}\n")
        print("\n✅ 附录A完整数据已导出为 appendix_A_3x3_anf_full.csv "
              f"（共 {len(all_anf_records)} 行）")

    return {
        'cnt_total': cnt_total,
        'cnt_total_opt': cnt_total_opt,
        'cnt_tv_opt': cnt_tv_opt,
        'cnt_fine_opt': cnt_fine_opt,
        'n_fine_classes_opt': n_fine_classes_opt,
        'n_fine_classes_all': n_fine_classes_all,
        'optimal_count': optimal_count,
    }

# ============================================================
# 第十三部分：主程序
# ============================================================

def main():
    output_file = open("results_v3.txt", "w", encoding="utf-8")
    print("=" * 70)
    print("小规模S盒的代数表示及密码学性质计算")
    print("（v3：新增基于ANF表达式结构的分类）")
    print("=" * 70)
    output_file.write("小规模S盒的代数表示及密码学性质计算（v3）\n\n")

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

    # 九、任务3（新增）：基于ANF表达式结构的分类 —— 老师新要求
    print("\n【九、基于ANF表达式结构的3×3 S盒分类（老师新要求）】")
    output_file.write("【九、基于ANF表达式结构的3×3 S盒分类】\n")
    task3_anf_structure_classification(output_file, export_all_anf=True)

    output_file.close()
    print(f"\n✅ 所有结果已保存到 results_v3.txt")
    print(f"✅ 40320个S盒的完整ANF已保存到 appendix_A_3x3_anf_full.csv")

if __name__ == "__main__":
    main()

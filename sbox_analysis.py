

from itertools import permutations
import time

# ============================================================
# 第一部分：基础工具函数
# ============================================================

def int_to_bits(val, n):
    """将整数转换为n位二进制列表，高位在前"""
    return [(val >> (n - 1 - i)) & 1 for i in range(n)]

def bits_to_int(bits):
    """将二进制列表转换为整数"""
    val = 0
    for b in bits:
        val = (val << 1) | b
    return val

def hamming_weight(val):
    """计算整数的汉明重量（二进制中1的个数）"""
    count = 0
    while val:
        count += val & 1
        val >>= 1
    return count

def inner_product(a, b, n):
    """计算两个n比特整数的内积（模2）"""
    return hamming_weight(a & b) % 2

# ============================================================
# 第二部分：S盒拆分与ANF求解
# ============================================================

def sbox_to_component_truth_tables(sbox, n, m):
    """
    将S盒拆分为各分量布尔函数的真值表
    输入：sbox = S盒真值表列表，n = 输入位数，m = 输出位数
    输出：长度为m的列表，每个元素是一个长度为2^n的真值表
    """
    size = 1 << n
    components = []
    for bit in range(m):
        tt = []
        for x in range(size):
            # 提取输出的第bit位（从高位开始，bit=0是最高位）
            y = sbox[x]
            bit_val = (y >> (m - 1 - bit)) & 1
            tt.append(bit_val)
        components.append(tt)
    return components

def mobius_transform(truth_table, n):
    """
    Möbius变换：从真值表计算ANF系数
    输入：truth_table = 长度为2^n的真值表，n = 变量个数
    输出：长度为2^n的ANF系数列表
    """
    size = 1 << n
    # 复制真值表
    anf = truth_table[:]
    # 蝶形运算（类似快速Walsh变换，但在GF(2)上）
    for i in range(n):
        step = 1 << i
        for j in range(size):
            if j & step:
                anf[j] ^= anf[j ^ step]
    return anf

def anf_to_string(anf_coeffs, n, var_names=None):
    """
    将ANF系数转换为可读的代数表达式字符串
    """
    if var_names is None:
        var_names = [f"x{i}" for i in range(n)]
    
    terms = []
    size = 1 << n
    for k in range(size):
        if anf_coeffs[k] == 1:
            if k == 0:
                terms.append("1")
            else:
                factors = []
                for i in range(n):
                    if (k >> i) & 1:
                        factors.append(var_names[i])
                terms.append("".join(factors))
    
    if not terms:
        return "0"
    return " ⊕ ".join(terms)

def compute_algebraic_degree(anf_coeffs, n):
    """计算布尔函数的代数次数"""
    max_deg = 0
    size = 1 << n
    for k in range(size):
        if anf_coeffs[k] == 1:
            deg = hamming_weight(k)
            if deg > max_deg:
                max_deg = deg
    return max_deg

def count_anf_terms(anf_coeffs):
    """计算ANF中非零项的个数"""
    return sum(anf_coeffs)

# ============================================================
# 第三部分：Walsh-Hadamard变换与非线性度
# ============================================================

def walsh_hadamard_transform(truth_table, n):
    """
    计算布尔函数的Walsh-Hadamard变换谱
    输入：truth_table = 长度为2^n的真值表
    输出：长度为2^n的Walsh谱值列表
    """
    size = 1 << n
    # 将 {0,1} 真值表转换为 {1,-1}
    f = [1 - 2 * truth_table[x] for x in range(size)]
    
    # 快速Walsh-Hadamard变换（蝶形运算）
    wht = f[:]
    for i in range(n):
        step = 1 << i
        for j in range(size):
            if not (j & step):
                u = wht[j]
                v = wht[j | step]
                wht[j] = u + v
                wht[j | step] = u - v
    return wht

def nonlinearity_boolean(truth_table, n):
    """计算单个布尔函数的非线性度"""
    wht = walsh_hadamard_transform(truth_table, n)
    max_abs = max(abs(w) for w in wht)
    return (1 << (n - 1)) - max_abs // 2

def nonlinearity_sbox(sbox, n, m):
    """
    计算S盒的非线性度
    遍历所有非零输出掩码v，计算组合分量函数的非线性度，取最小值
    """
    size = 1 << n
    min_nl = float('inf')
    
    for v in range(1, 1 << m):
        # 计算组合分量函数 g_v(x) = v · S(x) 的真值表
        tt = []
        for x in range(size):
            # v · S(x) 的内积
            tt.append(inner_product(v, sbox[x], m))
        nl = nonlinearity_boolean(tt, n)
        if nl < min_nl:
            min_nl = nl
    
    return min_nl

# ============================================================
# 第四部分：差分均匀性
# ============================================================

def difference_distribution_table(sbox, n):
    """
    构造差分分布表（DDT）
    输出：2^n × 2^n 的二维列表
    """
    size = 1 << n
    ddt = [[0] * size for _ in range(size)]
    
    for dx in range(size):
        for x in range(size):
            dy = sbox[x] ^ sbox[x ^ dx]
            ddt[dx][dy] += 1
    
    return ddt

def differential_uniformity(sbox, n):
    """计算S盒的差分均匀性"""
    ddt = difference_distribution_table(sbox, n)
    size = 1 << n
    max_val = 0
    for dx in range(1, size):  # 排除dx=0
        for dy in range(size):
            if ddt[dx][dy] > max_val:
                max_val = ddt[dx][dy]
    return max_val

# ============================================================
# 第五部分：线性逼近表与线性逼近优势
# ============================================================

def linear_approximation_table(sbox, n, m):
    """
    构造线性逼近表（LAT）
    输出：2^n × 2^m 的二维列表
    """
    size_in = 1 << n
    size_out = 1 << m
    lat = [[0] * size_out for _ in range(size_in)]
    
    for a in range(size_in):
        for b in range(size_out):
            count = 0
            for x in range(size_in):
                if inner_product(a, x, n) == inner_product(b, sbox[x], m):
                    count += 1
            lat[a][b] = count - (size_in >> 1)
    
    return lat

def max_linear_approximation(sbox, n, m):
    """计算S盒的最大线性逼近优势"""
    lat = linear_approximation_table(sbox, n, m)
    size_in = 1 << n
    size_out = 1 << m
    max_val = 0
    for a in range(size_in):
        for b in range(1, size_out):  # 排除b=0
            if abs(lat[a][b]) > max_val:
                max_val = abs(lat[a][b])
    return max_val

# ============================================================
# 第六部分：雪崩特性与不动点
# ============================================================

def fixed_points(sbox, n):
    """计算S盒的不动点个数"""
    count = 0
    for x in range(1 << n):
        if sbox[x] == x:
            count += 1
    return count

def sac_distance(sbox, n, m):
    """
    计算S盒各分量布尔函数的严格雪崩距离
    输出：长度为m的列表
    """
    size = 1 << n
    components = sbox_to_component_truth_tables(sbox, n, m)
    distances = []
    
    for comp_idx in range(m):
        tt = components[comp_idx]
        max_dist = 0
        for bit in range(n):
            e = 1 << (n - 1 - bit)
            count = 0
            for x in range(size):
                if tt[x] ^ tt[x ^ e]:
                    count += 1
            dist = abs(count - (size >> 1))
            if dist > max_dist:
                max_dist = dist
        distances.append(max_dist)
    
    return distances

# ============================================================
# 第七部分：完整S盒分析函数
# ============================================================

def analyze_sbox(sbox, n, m, name="S-box"):
    """对一个S盒进行完整的密码学指标分析"""
    results = {}
    results['name'] = name
    results['n'] = n
    results['m'] = m
    results['truth_table'] = sbox
    
    # 1. 拆分为分量布尔函数
    components = sbox_to_component_truth_tables(sbox, n, m)
    
    # 2. 计算各分量的ANF
    anf_list = []
    anf_strings = []
    degrees = []
    term_counts = []
    var_names = [f"x{i}" for i in range(n)]
    
    for i in range(m):
        anf = mobius_transform(components[i], n)
        anf_list.append(anf)
        anf_str = anf_to_string(anf, n, var_names)
        anf_strings.append(anf_str)
        deg = compute_algebraic_degree(anf, n)
        degrees.append(deg)
        tc = count_anf_terms(anf)
        term_counts.append(tc)
    
    results['anf_strings'] = anf_strings
    results['degrees'] = degrees
    results['sbox_degree'] = max(degrees)
    results['term_counts'] = term_counts
    
    # 3. 非线性度
    nl = nonlinearity_sbox(sbox, n, m)
    results['nonlinearity'] = nl
    
    # 4. 差分均匀性
    du = differential_uniformity(sbox, n)
    results['diff_uniformity'] = du
    
    # 5. 最大线性逼近优势
    mla = max_linear_approximation(sbox, n, m)
    results['max_linear_approx'] = mla
    
    # 6. 不动点
    fp = fixed_points(sbox, n)
    results['fixed_points'] = fp
    
    # 7. 雪崩距离
    sac_dist = sac_distance(sbox, n, m)
    results['sac_distances'] = sac_dist
    
    return results

def print_results(results, file=None):
    """格式化输出分析结果"""
    lines = []
    lines.append("=" * 70)
    lines.append(f"S盒名称：{results['name']}")
    lines.append(f"规模：{results['n']}×{results['m']}")
    lines.append(f"真值表：{[hex(x) for x in results['truth_table']]}")
    lines.append("-" * 70)
    
    lines.append("【代数正规型（ANF）】")
    for i in range(results['m']):
        lines.append(f"  y{i} = {results['anf_strings'][i]}")
    
    lines.append(f"\n【代数次数】")
    for i in range(results['m']):
        lines.append(f"  deg(y{i}) = {results['degrees'][i]}，项数 = {results['term_counts'][i]}")
    lines.append(f"  S盒代数次数 = {results['sbox_degree']}")
    
    lines.append(f"\n【非线性度】 NL = {results['nonlinearity']}")
    lines.append(f"【差分均匀性】 δ = {results['diff_uniformity']}")
    lines.append(f"【最大线性逼近值】 |LAT_max| = {results['max_linear_approx']}")
    lines.append(f"【最大线性逼近优势】 A_max = {results['max_linear_approx']}/{1 << results['n']}")
    lines.append(f"【不动点个数】 {results['fixed_points']}")
    lines.append(f"【严格雪崩距离】 {results['sac_distances']}")
    lines.append("=" * 70)
    
    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n\n")

# ============================================================
# 第八部分：3×3双射S盒的批量枚举与分类
# ============================================================

def enumerate_3x3_bijective_sboxes():
    """
    枚举所有3×3双射S盒（即8个元素的全排列），计算密码学指标并分类
    共有 8! = 40320 个
    """
    print("\n" + "=" * 70)
    print("正在枚举所有3×3双射S盒（共40320个）...")
    print("=" * 70)
    
    start_time = time.time()
    
    # 按（非线性度，差分均匀性，代数次数）进行分类
    categories = {}
    total = 0
    
    for perm in permutations(range(8)):
        sbox = list(perm)
        total += 1
        
        # 计算三个核心指标
        nl = nonlinearity_sbox(sbox, 3, 3)
        du = differential_uniformity(sbox, 3)
        
        # 代数次数
        components = sbox_to_component_truth_tables(sbox, 3, 3)
        max_deg = 0
        for comp in components:
            anf = mobius_transform(comp, 3)
            deg = compute_algebraic_degree(anf, 3)
            if deg > max_deg:
                max_deg = deg
        
        key = (nl, du, max_deg)
        if key not in categories:
            categories[key] = {'count': 0, 'example': sbox}
        categories[key]['count'] += 1
        
        if total % 5000 == 0:
            print(f"  已处理 {total}/40320 ...")
    
    elapsed = time.time() - start_time
    print(f"\n枚举完成，耗时 {elapsed:.1f} 秒")
    
    return categories, total

def print_3x3_classification(categories, total, file=None):
    """输出3×3 S盒的分类结果"""
    lines = []
    lines.append("\n" + "=" * 70)
    lines.append("3×3双射S盒分类结果（按非线性度、差分均匀性、代数次数分类）")
    lines.append(f"总数：{total}")
    lines.append("=" * 70)
    lines.append(f"{'非线性度':>8} {'差分均匀性':>10} {'代数次数':>8} {'数量':>8} {'占比':>8} {'示例S盒'}")
    lines.append("-" * 70)
    
    for key in sorted(categories.keys(), reverse=True):
        nl, du, deg = key
        count = categories[key]['count']
        example = categories[key]['example']
        ratio = f"{count/total*100:.2f}%"
        lines.append(f"{nl:>8} {du:>10} {deg:>8} {count:>8} {ratio:>8} {example}")
    
    # 统计最优S盒
    optimal_keys = [k for k in categories if k[0] == 2 and k[1] == 2 and k[2] == 2]
    optimal_count = sum(categories[k]['count'] for k in optimal_keys)
    lines.append("-" * 70)
    lines.append(f"最优S盒（NL=2, δ=2, deg=2）数量：{optimal_count}，占比：{optimal_count/total*100:.2f}%")
    lines.append("=" * 70)
    
    text = "\n".join(lines)
    print(text)
    if file:
        file.write(text + "\n\n")

# ============================================================
# 第九部分：主程序
# ============================================================

def main():
    output_file = open("results.txt", "w", encoding="utf-8")
    
    print("=" * 70)
    print("小规模S盒的代数表示及密码学性质计算")
    print("=" * 70)
    output_file.write("小规模S盒的代数表示及密码学性质计算\n\n")
    
    # ---- 3×3 S盒分析 ----
    print("\n【一、3×3 S盒分析】\n")
    output_file.write("【一、3×3 S盒分析】\n\n")
    
    # S1
    s1 = [5, 3, 6, 1, 0, 7, 2, 4]
    r1 = analyze_sbox(s1, 3, 3, "S1 (论文实例1)")
    print_results(r1, output_file)
    
    # S2
    s2 = [0, 1, 3, 6, 7, 4, 5, 2]
    r2 = analyze_sbox(s2, 3, 3, "S2 (论文实例2)")
    print_results(r2, output_file)
    
    # ---- 4×4 S盒分析 ----
    print("\n【二、4×4 S盒分析】\n")
    output_file.write("【二、4×4 S盒分析】\n\n")
    
    # PRESENT S盒
    present = [0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD,
               0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2]
    r_present = analyze_sbox(present, 4, 4, "PRESENT S盒")
    print_results(r_present, output_file)
    
    # GIFT S盒
    gift = [0x1, 0xA, 0x4, 0xC, 0x6, 0xF, 0x3, 0x9,
            0x2, 0xD, 0xB, 0x7, 0x5, 0x0, 0x8, 0xE]
    r_gift = analyze_sbox(gift, 4, 4, "GIFT S盒")
    print_results(r_gift, output_file)
    
    # Midori Sb0
    midori = [0xC, 0xA, 0xD, 0x3, 0xE, 0xB, 0xF, 0x7,
              0x8, 0x9, 0x1, 0x5, 0x0, 0x2, 0x4, 0x6]
    r_midori = analyze_sbox(midori, 4, 4, "Midori Sb0")
    print_results(r_midori, output_file)
    
    # PRESENT S盒的差分分布表
    print("\n【三、PRESENT S盒差分分布表】\n")
    output_file.write("【三、PRESENT S盒差分分布表】\n\n")
    ddt = difference_distribution_table(present, 4)
    header = "Δx\\Δy " + " ".join(f"{i:>3X}" for i in range(16))
    print(header)
    output_file.write(header + "\n")
    for dx in range(16):
        row = f"  {dx:>3X}  " + " ".join(f"{ddt[dx][dy]:>3}" for dy in range(16))
        print(row)
        output_file.write(row + "\n")
    print()
    output_file.write("\n")
    
    # PRESENT S盒的线性逼近表
    print("\n【四、PRESENT S盒线性逼近表】\n")
    output_file.write("【四、PRESENT S盒线性逼近表】\n\n")
    lat = linear_approximation_table(present, 4, 4)
    header = " a\\b  " + " ".join(f"{i:>3X}" for i in range(16))
    print(header)
    output_file.write(header + "\n")
    for a in range(16):
        row = f"  {a:>3X}  " + " ".join(f"{lat[a][b]:>3}" for b in range(16))
        print(row)
        output_file.write(row + "\n")
    print()
    output_file.write("\n")
    
    # ---- 4×4 S盒综合对比表 ----
    print("\n【五、4×4 S盒综合指标对比】\n")
    output_file.write("【五、4×4 S盒综合指标对比】\n\n")
    
    comparison = [r_present, r_gift, r_midori]
    header = f"{'S盒':>12} {'非线性度':>8} {'差分均匀性':>10} {'代数次数':>8} {'|LAT_max|':>10} {'不动点':>6}"
    print(header)
    output_file.write(header + "\n")
    print("-" * 60)
    output_file.write("-" * 60 + "\n")
    for r in comparison:
        row = f"{r['name']:>12} {r['nonlinearity']:>8} {r['diff_uniformity']:>10} {r['sbox_degree']:>8} {r['max_linear_approx']:>10} {r['fixed_points']:>6}"
        print(row)
        output_file.write(row + "\n")
    print()
    output_file.write("\n")
    
    # ---- 3×3 S盒批量枚举 ----
    print("\n【六、3×3双射S盒批量枚举与分类】")
    output_file.write("【六、3×3双射S盒批量枚举与分类】\n")
    
    categories, total = enumerate_3x3_bijective_sboxes()
    print_3x3_classification(categories, total, output_file)
    
    output_file.close()
    print(f"\n所有结果已保存到 results.txt 文件中。")

if __name__ == "__main__":
    main()
